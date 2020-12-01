import argparse
import os
import logging

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("file_list", help='File containing the Haddock-scores for each PDB')
    parser.add_argument("input_pdblist", help='Filtered PDB List containing the full path of the PDB')
    parser.add_argument("--top", help='After ranking, how many models should be considered when counting contacts, default=100', type=int, default=100)
    args = parser.parse_args()

    logging.basicConfig(level='DEBUG',
                        format='%(asctime)s L%(lineno)d %(levelname)s - %(message)s',
                        datefmt='%d/%m/%Y %H:%M:%S')

    # Load the haddock-scores
    haddockscore_dic = {}
    with open(args.file_list, 'r') as file_h:
        for line in file_h.readlines():
            haddock_score = line.split()[-2]
            pdb_name = line.split()[0].split(':')[-1][:-1]
            haddockscore_dic[pdb_name] = float(haddock_score)

    if len(haddockscore_dic) < args.top:
        logging.warning('Your input contains less than %i models, all models will be used.', args.top)

    top_models = sorted(haddockscore_dic)[:args.top]

    # Load the PDBs, note that here we are using it0 models
    logging.info('Extracting contacts from the top %i models', len(top_models))
    contacts = {}
    with open(args.input_pdblist, 'r') as fh:
        for pdb in fh.readlines():
            pdb = pdb[:-1]  # remove /n
            pdb_name = pdb.split('/')[-1]
            if pdb_name not in top_models:
                continue

            contact_file = pdb.replace('.pdb', '.contacts')

            if not os.path.isfile(contact_file):
                logging.warning('Contact file for %s not found! expected: %s', pdb, contact_file)

            with open(contact_file, 'r') as con_fh:
                for line in con_fh.readlines():
                    data = line.split()
                    chain_i = str(data[1])
                    chain_j = str(data[4])
                    resnum_i = int(data[0])
                    resnum_j = int(data[3])
                    pair = resnum_i, resnum_j

                    if chain_i not in contacts:
                        contacts[chain_i] = []
                    if chain_j not in contacts:
                        contacts[chain_j] = []

                    contacts[chain_i].append(resnum_i)
                    contacts[chain_j].append(resnum_j)
    
    logging.info('Result of the contact analysis')
    for chain in contacts:
        contact_str = ','.join(map(str, set(contacts[chain])))
        logging.info('Chain %s - %s', chain, contact_str)

    # done
