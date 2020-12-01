# load the models and the template, then prepare the submission
import argparse
import pandas as pd
import itertools
import subprocess
import os
import copy
import logging
import shutil


def get_models(path):
    """ Find the selection inside the scoring folder structure

    :param pdb_l: List containing the selection as "01_01"
    :param path: Root location of the scoring
    :type pdb_l: list
    :type path: string
    :return: Actual filename of the structure, ex: target164-scoring_0921_conv.pdb
    :rtype: list
    """
    # relate models with their cluster
    file_dic = {}
    with open(f'{path}/structures/it1/water/file.list', 'r') as fh:
        for idx, element in enumerate(fh.readlines()):
            score = float(element.split()[-2])
            pdb = element.split()[0].split(':')[-1][:-1]
            file_dic[idx + 1] = (pdb, score)
            # file_list.append((idx, pdb, score))

    cluster_dic = {}
    with open(f'{path}/structures/it1/water/analysis/cluster.out') as fh:
        for line in fh.readlines():
            cluster_id = int(line.split()[1])
            cluster_elements = list(map(int, line.split()[4:]))
            cluster_dic[cluster_id] = cluster_elements

    # add cluster information to file_list
    #  clusters are ranked by its top4 elements
    cluster_scores = []
    for cluster_id in cluster_dic:
        scores = [file_dic[cluster_element][-1] for cluster_element in cluster_dic[cluster_id]][:4]
        mean_score = sum(scores) / len(scores)
        cluster_scores.append((cluster_id, mean_score))

    cluster_scores.sort(key=lambda x: x[1])

    structure_list = []
    for overall_cluster_ranking, cluster_id in enumerate(cluster_scores):
        cluster_id, _ = cluster_id
        overall_cluster_ranking += 1
        for internal_cluster_ranking, cluster_element in enumerate(cluster_dic[cluster_id]):
            internal_cluster_ranking += 1
            pdb = file_dic[cluster_element][0]
            single_structure_ranking = [e for e in file_dic if file_dic[e][0] == pdb][0]
            structure_list.append((pdb, single_structure_ranking, overall_cluster_ranking, internal_cluster_ranking))
    
    done_pdbs = [e[0] for e in structure_list]

    for single_structure_ranking in file_dic:
        pdb, _ = file_dic[single_structure_ranking]
        if pdb not in done_pdbs:
            structure_list.append((pdb, single_structure_ranking, float('nan'), float('nan')))
    
    return structure_list


def pad_line(line):
    """Helper function to pad line to 80 characters in case it is shorter, borrowed from from pdb-tools (:

    :param line: Line extracted from a PDB file
    :type line: string
    :return: Padded line
    :rtype: string
    """
    size_of_line = len(line)
    if size_of_line < 80:
        padding = 80 - size_of_line + 1
        line = line.strip('\n') + ' ' * padding + '\n'
    return line[:81]  # 80 + newline character


def identify_chains(pdb_f):
    """" Read PDB structure and return its chainIDs

    :param pdb_f: PDB filename
    :type pdb_f: string
    :return: List with chains found on the input PDB
    :rtype: list
    """
    chainid_l = []
    with open(pdb_f) as fh:
        for line in fh.readlines():
            if line.startswith('ATOM'):
                line = pad_line(line)
                chainid = line[21].strip()[:1]
                if not chainid.isdigit():  # ignore NUMERICAL chains
                    chainid_l.append(chainid)

    chainid_l = list(set(chainid_l))
    # WARNING: This expects chains to be sequential!
    # chain1 = A, chain2 = B or chain1 = D, chain = E, etc
    chainid_l.sort()
    return chainid_l


def load_seq(pdb_f):
    """ Load a PDB filename into a dictionary {chain: {resnum: resname}}

    :param pdb_f: PDB file name
    :type pdb_f: string
    :return: {chain: {resnum: resname}}
    :rtype: dict
    """
    aa_dic = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'DC': 'C', 'DA': 'A', 'DG': 'G', 'DT': 'T',
              'ADE': 'A', 'THY': 'T', 'GUA': 'G', 'CYT': 'C'}
    seq_dic = {}
    for l in open(pdb_f):
        if 'ATOM' in l[:4]:
            chain = l[21]
            resnum = int(l[22:26])
            resname = l[17:20].split()[0]
            if chain not in seq_dic:
                seq_dic[chain] = {}
            try:
                name = aa_dic[resname]
            except KeyError:
                name = 'X'
            seq_dic[chain][resnum] = name

    return seq_dic


def replace_chain(pdbf, old_chain, new_chain, overwrite=True):
    """ Replace a specified chain with a new one in a PDB file

    :param pdbf: PDB filename
    :param old_chain: Old chain
    :param new_chain: New chain
    :param overwrite: Overwrite the PDB or create a new one
    :type pdbf: string
    :type old_chain: string
    :type new_chain: string
    :type overwrite: bool
    :return: Replaced PDB filename
    :rtype: string
    """
    new_pdb = pdbf + '_'
    with open(new_pdb, 'w') as out_fh:
        with open(pdbf) as in_fh:
            for line in in_fh.readlines():
                if line.startswith('ATOM'):
                    current_chain = line[21]
                    if current_chain == old_chain:
                        line = line[:21] + new_chain + line[22:]
                out_fh.write(line)
        in_fh.close()
    out_fh.close()

    if overwrite:
        os.rename(new_pdb, pdbf)
        return pdbf
    else:
        return new_pdb


def renumber(pdbf, renumber_dic, target_chain, overwrite=True):
    """ Renumber a specific chain of a PDB file based on a reference dictionary

    :param pdbf: PDB filename
    :param renumber_dic: Reference dicionary {template_resnum: target_resnum}
    :param target_chain: Chain to be renumbered
    :param overwrite: Overwrite the PDB or create a new one
    :type pdbf: string
    :type renumber_dic: dict
    :type target_chain: string
    :type overwrite: bool
    :return: The renumbered PDB filename
    :rtype: string
    """
    new_pdb = pdbf + '_'
    ignored_res = []
    with open(new_pdb, 'w') as out_fh:
        with open(pdbf) as in_fh:
            for line in in_fh.readlines():
                if line.startswith('ATOM'):
                    current_res = int(line[22:26])
                    current_chain = line[21]
                    if current_chain == target_chain:
                        try:
                            new_res = renumber_dic[current_res]
                            line = line[:22] + '{:>4}'.format(new_res) + line[26:]
                            out_fh.write(line)
                        except:
                            # Residue not found in reference, IGNORE
                            ignored_res.append(current_res)
                    else:
                        out_fh.write(line)
                else:
                    out_fh.write(line)
        in_fh.close()
    out_fh.close()

    if ignored_res:
        ignored_res_str = ', '.join(map(str, list(set(ignored_res))))
        logging.warning(f'{pdbf} Chain {target_chain} Res {ignored_res_str} not found in reference, discarded.')

    if overwrite:
        os.rename(new_pdb, pdbf)
        return pdbf
    else:
        return new_pdb


def match(pdb_l, template):
    """ Use sequence alignment to match the selection to the provided template by calculating the identity of them
        with a combinatorial product

    :param pdb_l: List containing multiple PDB filenames
    :param template: CAPRI .brk template filename
    :type pdb_l: list
    :type template: string
    :return: List of the matched/renumbered PDB filenames
    :rtype: list
    """
    # TODO: Improve the readability of this function
    match_l = []
    template_seq_dic = load_seq(template)
    template_chains = identify_chains(template)
    for pdb in pdb_l:
        pdb_seq_dic = load_seq(pdb)
        pdb_chains = identify_chains(pdb)
        identity_dic = {}
        for ref_chain, target_chain in itertools.product(template_chains, pdb_chains):

            ref_seq = ''.join(list(template_seq_dic[ref_chain].values()))
            target_seq = ''.join(list(pdb_seq_dic[target_chain].values()))

            # TODO: use python's tempfile
            open('seq.fasta', 'w').write(f'>ref\n{ref_seq}\n>target\n{target_seq}\n')

            cmd = '/Users/rodrigo/software/clustalo -i seq.fasta --outfmt=clu --resno --wrap=9000 --force'
            # logging.debug(f'Running align command {cmd}')
            p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out = p.communicate()
            os.remove('seq.fasta')

            aln_data = out[0].decode('utf-8').split()
            ref_aln = aln_data[6]
            target_aln = aln_data[9]
            counter_a = 0
            counter_b = 0
            numbering_dic = {}
            for i in range(len(ref_aln)):
                ref_char = ref_aln[i]
                target_char = target_aln[i]
                try:
                    ref_resnum = list(template_seq_dic[ref_chain])[counter_a]
                except IndexError:
                    # Ref sequence exhausted, ignore
                    ref_resnum = '-'
                try:
                    target_resnum = list(pdb_seq_dic[target_chain])[counter_b]
                except IndexError:
                    # Target sequence exhausted, ignore
                    target_resnum = '-'

                if '-' not in ref_char:
                    counter_a += 1
                if '-' not in target_char:
                    counter_b += 1
                if '-' not in ref_char and '-' not in target_char:
                    numbering_dic[target_resnum] = ref_resnum

            identity = out[0].decode('utf-8').count('*') / float(len(ref_seq))
            coverage = len(numbering_dic) / len(ref_aln)

            logging.debug(f"{ref_chain}, {target_chain}, {identity}, {coverage}")
            # logging.debug(f'>R:{ref_chain}\n{ref_aln}')
            # logging.debug(f'>T:{target_chain}\n{target_aln}')

            try:
                identity_dic[ref_chain].append((target_chain, identity, coverage, numbering_dic))
            except KeyError:
                identity_dic[ref_chain] = [(target_chain, identity, coverage, numbering_dic)]

        for i, ref_c in enumerate(template_chains):
            target_info_list = [(v[0], v[1], v[2]) for v in identity_dic[ref_c]]
            # sort by identity and coverage
            sorted_target_list = sorted(target_info_list, key=lambda x: (-x[2], x[1]))
            # create a catalog with possible numbering references
            numbering_dic_catalog = dict([(v[0], v[3]) for v in identity_dic[ref_c]])

            if len(set([e[1] for e in sorted_target_list])) == 1:
                # this is a homo-something, match is sequentialy
                selected_chain = pdb_chains[i]
            else:
                # get the highest identity/coverage
                selected_chain = sorted_target_list[0][0]

            # select the correct numbering dictionary
            selected_numbering_dic = numbering_dic_catalog[selected_chain]
            # just for readability:
            old_chain = selected_chain
            new_chain = ref_c

            # replace the target chain (old) with the same observed in the reference (new)
            chain_matched_pdb = replace_chain(pdb, old_chain, new_chain, overwrite=True)

            # renumber!
            #  Note, if residue is present in target and not in reference it will be DELETED, use with caution
            renumbered_pdb = renumber(chain_matched_pdb, selected_numbering_dic, new_chain, overwrite=True)

            if renumbered_pdb not in match_l:
                match_l.append(renumbered_pdb)

    return match_l


def create_ensemble(pdb_l, ensemble_name):
    """ Join a list of PDBs to create an ensemble with pdb_mkensemble

    :param pdb_l: List containing the PDB filenames to be made into an ensemble
    :param ensemble_name: Filename of the ensemble to be generated
    :type pdb_l: list 
    :type ensemble_name: string  
    :return: True if ensemble was generated
    :rtype: bool
    """
    pdb_str = ' '.join(pdb_l)
    cmd = f'pdb_mkensemble {pdb_str}'
    logging.debug(f'Create ensemble with command {cmd}')
    output_f = open(ensemble_name, 'w')
    p = subprocess.Popen(cmd.split(), stdout=output_f, stderr=subprocess.PIPE)
    p.communicate()
    if not os.path.isfile(ensemble_name):
        logging.error(f'Could not generate the ensemble {ensemble_name}')
    return True


def add_csb_header(submission):
    """ Add the CSB Header into a string containg already the Capri Header and the PDB Ensemble

    :param submission: String containing the HEADER and the ENSEMBLE
    :type submission: string
    :return: Capri HEADER + CSB HEADER + Ensemble
    :rtype: string
    """
    # The CSB header must be added right after the COMPND remark
    pos = 0
    submission_aslist = submission.split('\n')
    for l in submission_aslist:
        if not l.startswith('COMPND'):
            continue
        else:
            pos += 1
    return '\n'.join(submission_aslist[:pos + 1] + [csb_header()] + submission_aslist[pos + 1:])


def csb_header():
    """ The CSB Header used for CAPRI Scorer Submission

    :return: The csb header
    :rtype: string
    """
    header = '''AUTHOR   1 Alexandre M.J.J. Bonvin, Zuzana Jandova, Francesco Ambrosetti,
AUTHOR   2 Jorge Roel Touris, Panagiotis Koukos, Charlotte van Noort, Siri van Keulen,
AUTHOR   3 Rodrigo Honorato, Brian Jimenez Garcia, Adrien Melquiond, Manon Reau
REMARK   1
REMARK   1 Models obtained with HADDOCK2.4
REMARK   2
REMARK   2 REFERENCE 1
REMARK   2  AUTH   S.J.de Vries, M. van Dijk, A.M.J.J. Bonvin
REMARK   2  TITL   The HADDOCK web server for data-driven biomolecular docking
REMARK   2  TITL 2 of HADDOCK2.0 on the CAPRI targets
REMARK   2  REF    Nature Protocols              V.5      883 2010
REMARK   3 REFERENCE 2
REMARK   3  AUTH   G.C.P van Zundert, J.P.G.L.M. Rodrigues, M. Trellet,
REMARK   3  AUTH 2 C. Schmitz, P.L. Kastritis, E. Karaca, A.S.J. Melquiond
REMARK   3  AUTH 3 M. van Dijk, S.J. de Vries, A.M.J.J. Bonvin
REMARK   3  TITL   The HADDOCK2.2 webserver: User-friendly integrative
REMARK   3  TITL 2 modeling of biomolecular complexes.
REMARK   3  REF    J.Mol.Biol., 428, 720-725 (2015).'''
    return ''.join(header)


def tidy(pdb_f):
    cmd = f'pdb_tidy {pdb_f}'
    logging.debug(f'Tidying up with command {cmd}')
    tidy_output_f = open(f'{pdb_f}_', 'w')
    p = subprocess.Popen(cmd.split(), stdout=tidy_output_f, stderr=subprocess.PIPE)
    p.communicate()
    if not os.path.isfile(f'{pdb_f}_'):
        logging.error(f'Could not generate the tidy PDB  {pdb_f}_')
    else:
        shutil.move(f'{pdb_f}_', pdb_f)
        return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('run_path', help="Location of the run")
    parser.add_argument('template', help="Template provided by CAPRI")
    args = parser.parse_args()

    if not os.path.isdir('selection'):
        os.mkdir('selection')

    run_path = args.run_path
    # run_path = '/Users/rodrigo/repos/capri-r51-target183/runs/31349-nsp8-hexadeca-exosome'

    template = args.template
    # template = '/Users/rodrigo/repos/capri-r51-target183/capri_51_183.brk'

    data = get_models(run_path)

    # make this into a dataframe
    df = pd.DataFrame(data, columns=['pdb', 'single_structure_ranking', 'overall_cluster_ranking', 'internal_cluster_ranking'])

    # now select top10 cluster then fill until we reach 100
    top1_top10_cluster = df[(df['overall_cluster_ranking'] <= 5) & (df['internal_cluster_ranking'] == 1)]['pdb']
    top2_top10_cluster = df[(df['overall_cluster_ranking'] <= 5) & (df['internal_cluster_ranking'] == 2)]['pdb']

    top10 = list(itertools.chain(*zip(top1_top10_cluster, top2_top10_cluster)))

    # sort the dataframe by single-structure ranking and fill until we get 100 models
    single_structure_sorted_df = df.sort_values(by='single_structure_ranking')

    selection = copy.deepcopy(top10)

    while len(selection) != 100:
        for pdb in single_structure_sorted_df['pdb']:
            if pdb not in selection:
                selection.append(pdb)
                break
    
    with open('selection/selection.txt', 'w') as fh:
        for pdb in selection:
            fh.write(f'{run_path}/structures/it1/water/{pdb}\n')
    fh.close()

    # =====================================================================================
    #
    # WARNING: This here is a hardcoded bypass specifically for CAPRI51 Target183
    #
    # =====================================================================================
    renumbered_pdb_list = []
    for pdb in selection:
        output_pdb = f'selection/{pdb}'
        intput_pdb = f'{run_path}/structures/it1/water/{pdb}'
        with open(output_pdb, 'w') as out_fh:
            with open(intput_pdb, 'r') as fh:
                for line in fh.readlines():
                    if line.startswith('ATOM'):
                        chain = line[21]
                        resnum = int(line[22:26])
                        if chain == 'A':
                            if resnum not in range(2282, 2397):
                                line = line[:21] + '0' + line[22:]
                            else:
                                line = line[:21] + 'X' + line[22:]
                        if chain == 'B':
                            if resnum not in range(576, 845):
                                line = line[:21] + '1' + line[22:]
                            else:
                                line = line[:21] + 'Y' + line[22:]
                        out_fh.write(line)
        out_fh.close()
        renumbered_pdb_list.append(output_pdb)
    # =====================================================================================
                
    matched_pdbs = match(renumbered_pdb_list, template)

    # logging.info(f'Creating HEADER')
    header_str = ''
    permitted_records = ('HEADER', 'COMPND', 'SEQRES', 'REMARK')
    with open(template) as fh:
        for l in fh.readlines():
            if l.startswith(permitted_records):
                header_str += l
            if l.startswith('MODEL'):
                break

    # logging.info(f'Creating the ensemble with pdb_mkensemble')
    ensemble_f = create_ensemble(matched_pdbs, 'temp_ens.pdb')

    logging.info('Merging HEADER with the ensemble')
    submission_str = header_str
    with open('temp_ens.pdb') as fh:
        for l in fh.readlines():
            if l.startswith('REMARK'):
                continue
            submission_str += l

    os.remove('temp_ens.pdb')

    logging.info('Adding the CSB Header')
    submission_str = add_csb_header(submission_str)
    with open('selection/Target183_selection.pdb', 'w') as fh:
        fh.write(submission_str)
    fh.close()

    logging.info('Adding TER with pdb_tidy')
    finalized = tidy('selection/Target183_selection.pdb')

    # done
