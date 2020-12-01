# Capri Round 51 (COVID19 thematic) Target 183

Selected models for this target are available at `/selection`


## Requirements
* Python 3x + Pandas
* [haddock-tools](https://github.com/haddocking/haddock-tools)
* [pdb-tools](https://github.com/haddocking/pdb-tools)
* [HADDOCK v2.4 Webserver File Interface](https://bianca.science.uu.nl/haddock2.4/submit_file) to replicate the results

***
## Information

```
Capri target description:
EXOSC8 is a non-catalytic component of the RNA exosome complex (9 subunits). PDB:2NN6 seems to contain the whole exosome. Nsp8-EXOSC2, Nsp8-EXOSC3 and Nsp8-EXOSC5 are also part of the list of Virus-Host interactions with 80 % template coverage. Predictors may need to consider interactions with these subunits as well. Note that Q13868 (EXOSC2) and Q9NPD3 (EXOS4) do not show up as interactors of Nsp7 and Nsp12 (which form a complex with NSsp18) in Gordon et al. dataset.
```

This could be interpreted as:
```
Nsp8 -EXOSC2 ✅
Nsp8 -EXOSC3 ✅
Nsp8 -EXOSC5 ✅
Nsp7 -EXOSC2 🚫
Nsp7 -EXOS4  🚫
```
***

## Approach

Dock the Nsp7/8 dimer against the full Exosome and then filter out conformations in which Nsp7 interacts with EXOSC2 and EXOSC4. Use the remaining conformations for a contact analysis and take the result of this contact analysis as restraints for another round of docking now using the hexadecamer of Nsp7/8 which was part of T184.

***

## Input

1. Mol A - Nsp7/Nsp8
    The PDB [6m5i](https://www.rcsb.org/structure/6M5I) is Nsp7/Nsp8 complex but the structure has a very bad quality and it’s also incomplete, it's a perfect match with Nsp8. We run molecular refinement to improve its quality - [Nsp7/Nsp8 refinement](runs/nsp7_8-refinement.hson)

2. Mol B - Exosome

    The PDB [2nn6](https://www.rcsb.org/structure/2NN6) has the whole exosome but the structure is also not of good quality. Refine it as well - [Exosome refinement](runs/exosome-refinement.json)


## Processing

To use it in HADDOCK the Top1 of each of the refinement run was renumbered/rechained using [pdb-tools](http://github.com/haddocking/pdb-tools) and its solvent accessible surface calculated with [haddock-tools/calc-accessibility.py](http://github.com/haddocking/haddock-tools)

```
# Renumbered residue mapping
NSP8 A  1-115
NSP7 B  116-196

EXOSC9  A
EXOSC4  B   304-537
EXOSC8  C   538-807
EXOSC5  D   808-1012
EXOSC7  E
EXOSC6  F
EXOSC3  G   1511-1746
EXOSC2  H   1747-1996
EXOSC1  I
```
```
$ python haddock-tools/calc-accessibility.py input/nsp8_A-nsp7_B.pdb
$ python haddock-tools/calc-accessibility.py input/exosome.pdb
```

### Get Accessible residues of Nsp8

Use some simple python scripting to find the intersection between all accessible residues in the dimer and those ones belonging only to Nsp8.


```python
$ python
>>> nsp8 = list(range(1,115))
>>> all_accessible = [1,2,3,4,5,6,7,9,10,13,14,18,20,21,23,25,28,29,32,35,36,37,42,47,48,49,51,58,60,63,64,67,68,69,70,72,75,77,79,81,82,86,87,88,89,90,92,93,94,95,97,98,99,100,102,103,105,107,113,115,116,117,119,120,123,126,130,133,136,138,139,140,141,142,145,146,148,149,152,153,156,158,159,160,161,162,165,178,179,182,184,185,188,189,191,192,193,194,195,196]
>>> set(nsp8).intersection(all_accessible)
{1, 2, 3, 4, 5, 6, 7, 9, 10, 13, 14, 18, 20, 21, 23, 25, 28, 29, 32, 35, 36, 37, 42, 47, 48, 49, 51, 58, 60, 63, 64, 67, 68, 69, 70, 72, 75, 77, 79, 81, 82, 86, 87, 88, 89, 90, 92, 93, 94, 95, 97, 98, 99, 100, 102, 103, 105, 107, 113}
```


### Get Accessible residues of EXOSC2/3/4/8

Same for EXOSC2/3/4/8:
```python
$ python
>>> exosc2 = list(range(1747,1996))
>>> exosc3 = list(range(1511,1746))
>>> exosc5 = list(range(808,1012))
>>> exosc8 = list(range(538,807))
>>> all_accessible = [3,8,11,12,15,16,19,20,21,25,26,28,31,34,36,42,43,52,60,62,64,65,67,68,69,70,72,75,77,89,92,93,94,96,97,98,99,111,114,117,118,119,120,123,125,126,127,164,166,168,169,170,171,173,175,176,177,178,180,181,182,184,188,200,201,203,218,227,228,229,238,241,243,244,245,249,252,253,256,259,260,263,264,266,267,270,271,273,274,275,277,278,279,280,281,282,286,287,290,291,292,293,294,297,299,301,303,305,306,307,310,311,312,314,316,317,320,322,324,329,341,356,357,358,359,360,362,363,364,365,379,380,381,383,384,386,387,388,389,390,391,393,405,410,411,413,414,427,450,463,464,471,472,475,476,477,478,480,489,500,502,504,505,507,508,511,512,515,518,519,522,523,525,526,527,530,533,534,536,537,538,539,540,541,542,543,544,547,550,551,552,553,556,557,559,560,562,565,566,567,569,574,577,586,597,598,599,601,602,603,605,606,607,609,611,623,625,627,628,629,645,648,651,653,654,659,660,661,693,695,698,700,702,703,704,705,706,707,709,710,711,713,714,715,716,717,719,721,722,733,744,747,748,750,759,760,762,775,777,778,781,784,785,787,788,791,792,795,798,799,801,802,803,805,806,807,808,810,819,822,834,847,848,849,850,852,854,861,863,865,866,867,871,877,880,881,887,889,902,907,926,938,939,940,941,942,947,948,949,952,953,955,957,965,966,975,977,981,984,985,988,989,991,992,995,999,1002,1006,1007,1010,1011,1013,1014,1016,1019,1020,1024,1029,1036,1039,1042,1044,1046,1048,1049,1051,1062,1063,1073,1075,1076,1078,1079,1080,1081,1083,1088,1090,1092,1099,1100,1101,1102,1103,1104,1105,1107,1108,1112,1118,1119,1122,1125,1129,1130,1136,1137,1138,1139,1172,1175,1177,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1196,1197,1198,1209,1210,1218,1219,1222,1225,1235,1237,1245,1246,1252,1255,1256,1259,1263,1265,1266,1267,1268,1270,1271,1274,1277,1278,1279,1281,1282,1284,1285,1286,1287,1288,1289,1290,1291,1292,1293,1297,1299,1301,1306,1309,1317,1318,1327,1329,1330,1331,1332,1333,1335,1337,1339,1348,1349,1351,1352,1353,1354,1356,1357,1358,1360,1375,1377,1378,1380,1381,1398,1429,1430,1431,1432,1433,1434,1435,1442,1444,1445,1446,1450,1460,1470,1471,1472,1473,1475,1476,1479,1482,1485,1486,1489,1493,1494,1497,1501,1504,1505,1506,1507,1508,1509,1510,1511,1512,1514,1515,1516,1517,1519,1520,1527,1529,1531,1532,1533,1534,1536,1539,1542,1544,1546,1556,1558,1560,1561,1563,1564,1565,1566,1567,1568,1570,1578,1583,1584,1585,1593,1594,1596,1597,1598,1601,1605,1615,1618,1619,1620,1621,1622,1623,1625,1626,1627,1629,1630,1631,1637,1639,1643,1644,1654,1655,1656,1657,1658,1659,1660,1662,1665,1666,1667,1672,1674,1680,1684,1685,1687,1689,1690,1693,1694,1696,1697,1699,1713,1714,1716,1717,1720,1723,1726,1727,1730,1732,1733,1734,1735,1737,1738,1742,1745,1746,1747,1753,1754,1756,1758,1759,1760,1762,1763,1765,1769,1770,1771,1772,1773,1783,1784,1785,1786,1787,1788,1792,1799,1801,1802,1803,1812,1814,1815,1816,1817,1821,1823,1825,1831,1832,1833,1836,1838,1839,1842,1843,1846,1855,1856,1858,1859,1860,1861,1864,1868,1869,1870,1871,1882,1886,1889,1890,1891,1892,1895,1897,1899,1901,1917,1919,1920,1921,1922,1923,1924,1925,1926,1927,1928,1929,1930,1931,1932,1949,1951,1953,1956,1960,1964,1967,1968,1970,1971,1972,1973,1975,1976,1978,1981,1982,1985,1988,1989,1990,1992,1993,1994,1995,1996,1997,2003,2004,2005,2008,2009,2010,2012,2014,2015,2017,2021,2022,2023,2033,2034,2035,2036,2037,2038,2039,2040,2041,2042,2047,2048,2049,2050,2051,2053,2054,2055,2056,2057,2058,2059,2060,2061,2063,2065,2068,2072,2073,2075,2080,2082,2083,2084,2086,2089,2090,2091,2093,2097,2099,2100,2104,2105,2106,2111,2114,2116,2122,2125,2128,2129,2130,2131,2139,2140,2141,2145,2146,2148,2150,2152,2154,2156,2157,2158,2159,2160,2162,2166,2167,2169,2171,2173,2176]
>>> selected_res = []
>>> selected_res += list(set(exosc2).intersection(all_accessible))
>>> selected_res += list(set(exosc3).intersection(all_accessible))
>>> selected_res += list(set(exosc5).intersection(all_accessible))
>>> selected_res += list(set(exosc8).intersection(all_accessible))
>>> selected_res
[1747, 1753, 1754, 1756, 1758, 1759, 1760, 1762, 1763, 1765, 1769, 1770, 1771, 1772, 1773, 1783, 1784, 1785, 1786, 1787, 1788, 1792, 1799, 1801, 1802, 1803, 1812, 1814, 1815, 1816, 1817, 1821, 1823, 1825, 1831, 1832, 1833, 1836, 1838, 1839, 1842, 1843, 1846, 1855, 1856, 1858, 1859, 1860, 1861, 1864, 1868, 1869, 1870, 1871, 1882, 1886, 1889, 1890, 1891, 1892, 1895, 1897, 1899, 1901, 1917, 1919, 1920, 1921, 1922, 1923, 1924, 1925, 1926, 1927, 1928, 1929, 1930, 1931, 1932, 1949, 1951, 1953, 1956, 1960, 1964, 1967, 1968, 1970, 1971, 1972, 1973, 1975, 1976, 1978, 1981, 1982, 1985, 1988, 1989, 1990, 1992, 1993, 1994, 1995, 1536, 1539, 1542, 1544, 1546, 1556, 1558, 1560, 1561, 1563, 1564, 1565, 1566, 1567, 1568, 1570, 1578, 1583, 1584, 1585, 1593, 1594, 1596, 1597, 1598, 1601, 1605, 1615, 1618, 1619, 1620, 1621, 1622, 1623, 1625, 1626, 1627, 1629, 1630, 1631, 1637, 1639, 1643, 1644, 1654, 1655, 1656, 1657, 1658, 1659, 1660, 1662, 1665, 1666, 1667, 1672, 1674, 1680, 1684, 1685, 1687, 1689, 1690, 1693, 1694, 1696, 1697, 1699, 1713, 1714, 1716, 1717, 1720, 1723, 1726, 1727, 1730, 1732, 1733, 1734, 1735, 1737, 1738, 1742, 1745, 1511, 1512, 1514, 1515, 1516, 1517, 1519, 1520, 1527, 1529, 1531, 1532, 1533, 1534, 902, 907, 926, 808, 810, 938, 939, 940, 941, 942, 819, 947, 948, 822, 949, 952, 953, 955, 957, 834, 965, 966, 847, 848, 849, 850, 975, 852, 977, 854, 981, 984, 985, 988, 861, 989, 863, 991, 865, 866, 867, 992, 995, 871, 999, 1002, 877, 1006, 1007, 880, 881, 1010, 1011, 887, 889, 538, 539, 540, 541, 542, 543, 544, 547, 550, 551, 552, 553, 556, 557, 559, 560, 562, 565, 566, 567, 569, 574, 577, 586, 597, 598, 599, 601, 602, 603, 605, 606, 607, 609, 611, 623, 625, 627, 628, 629, 645, 648, 651, 653, 654, 659, 660, 661, 693, 695, 698, 700, 702, 703, 704, 705, 706, 707, 709, 710, 711, 713, 714, 715, 716, 717, 719, 721, 722, 733, 744, 747, 748, 750, 759, 760, 762, 775, 777, 778, 781, 784, 785, 787, 788, 791, 792, 795, 798, 799, 801, 802, 803, 805, 806]
```

## Docking stage 1

* MolA = Nsp8 - Solvent accessible active
* MolB = EXOSC2/3/5/8 - Surface passive

Using "Optimize for Bioinformatics prediction" option

* [Docking parameters](runs/28513-nsp8-surf-act-exosc2_3_5-passive_ncvpart.json)

_Submit it via HADDOCK's file interface and download/uncompress in `runs/`_

## Analysis stage 1

1. Filter out Nsp7+EXOSC2/EXOS4 contacts 

This script will first calculate the contacts of each single model present inside the `it0` directory of the simulation using [haddock-tools/contact-chainID](http://github.com/haddocking/haddock-tools), this will generate many `.contacts` file inside the directory. Then for each model a function will check how many contacts are part of the "forbidden contacts" (lower than 20%), those being nsp7+exosc2/4. A text file will be written containing the filtered models.

```
$ python scripts/filter-contacts.py -h
usage: filter-contacts.py [-h] [--np NP] [--cutoff CUTOFF] run_directory contact_exec

positional arguments:
  run_directory    file.nam generated by HADDOCK
  contact_exec     Compiled contact script

optional arguments:
  -h, --help       show this help message and exit
  --np NP          Number of processors to use
  --cutoff CUTOFF  Cutoff of forbidden contacts allowed in the PDB to be filtered, float between 0 and 1

$ python scripts/filter-contacts.py runs/28513-nsp8-surf-act-exosc2_3_5-passive_ncvpart ~/repos/haddock-tools/contact-chainID --np 8
```

2. Run contact analysis

We then rank the selected models by their HADDOCK-score and extract residues that were observed to be in contact for a specific subset. Here we are using the Top10 models since Haddock did a good job of enriching the desired contacts.


```
$ python scripts/contact_analysis.py -h
usage: contact_analysis.py [-h] [--top TOP] file_list input_pdblist

positional arguments:
  file_list      File containing the Haddock-scores for each PDB
  input_pdblist  Filtered PDB List containing the full path of the PDB

optional arguments:
  -h, --help     show this help message and exit
  --top TOP      After ranking, how many models should be considered when counting contacts, default=100
$ python scripts/contact_analysis.py runs/28513-nsp8-surf-act-exosc2_3_5-passive_ncvpart/structures/it0/file.list filtered-pdbs.list
```

The contact analysis shows that basically all of Nsp8 makes contact with the Exosome, however some of the Exosome residues have been filtered out, from the initial 258 to 99.

```python
$ python
>>> filt_exosome = [1527,1545,1546,1547,1548,1556,1557,1558,1559,1560,1561,1562,1563,1564,1565,1566,1567,1568,547,548,1574,551,552,553,559,2105,1618,1619,1620,1621,1622,1623,1626,1627,606,2145,2146,1637,2152,2153,2154,2156,2158,2160,1653,2166,1655,1656,1657,1658,1659,2167,1660,1661,1654,2169,2171,533,1666,651,653,654,1680,1682,1683,1684,1685,1686,1687,1688,1689,1690,1691,1693,1694,1696,1697,1698,1699,1700,175,1713,1729,1730,1731,1732,1733,1734,1735,1736,1737,1738,1741,1742,1745,721,1747,1756,1759,1760,1761,1762,1763,1765,1769,1772,1782,1783,1784,1785,1786,1787,760,775,776,1801,1802,1803,777,778,781,780,784,290,292,293,295,296,297,298,299,300,301,302,1853,1854,1855,1856,1857,1858,1859,1868,1869,1870,1871,1872,1889,1890,1892,1895,1896,1897,1899,1919,1920,410,411,999,938,939,1002,1003,953,955,1994,1995,992,993,995,996,1511,1512,1513,1514,1515,1516,1517,1518,1519,1520,1006,1010,1007,1009,1012,1011,1521,1528,1529,1530,1531,1532,1533,1534,1526]
>>> exosc8 = list(range(538,807))
>>> sele_exosc8 = list(set(exosc8).intersection(filt_exosome))
>>> sele_exosc8.sort()
>>> sele_exosc8
[547, 548, 551, 552, 553, 559, 606, 651, 653, 654, 721, 760, 775, 776, 777, 778, 780, 781, 784]
```

Visualize this in PyMol with:
```
$ pymol input/exosome.pdb
PyMol> color white
PyMol> sele sas, resid 1747+1753+1754+1756+1758+1759+1760+1762+1763+1765+1769+1770+1771+1772+1773+1783+1784+1785+1786+1787+1788+1792+1799+1801+1802+1803+1812+1814+1815+1816+1817+1821+1823+1825+1831+1832+1833+1836+1838+1839+1842+1843+1846+1855+1856+1858+1859+1860+1861+1864+1868+1869+1870+1871+1882+1886+1889+1890+1891+1892+1895+1897+1899+1901+1917+1919+1920+1921+1922+1923+1924+1925+1926+1927+1928+1929+1930+1931+1932+1949+1951+1953+1956+1960+1964+1967+1968+1970+1971+1972+1973+1975+1976+1978+1981+1982+1985+1988+1989+1990+1992+1993+1994+1995+1536+1539+1542+1544+1546+1556+1558+1560+1561+1563+1564+1565+1566+1567+1568+1570+1578+1583+1584+1585+1593+1594+1596+1597+1598+1601+1605+1615+1618+1619+1620+1621+1622+1623+1625+1626+1627+1629+1630+1631+1637+1639+1643+1644+1654+1655+1656+1657+1658+1659+1660+1662+1665+1666+1667+1672+1674+1680+1684+1685+1687+1689+1690+1693+1694+1696+1697+1699+1713+1714+1716+1717+1720+1723+1726+1727+1730+1732+1733+1734+1735+1737+1738+1742+1745+1511+1512+1514+1515+1516+1517+1519+1520+1527+1529+1531+1532+1533+1534+902+907+926+808+810+938+939+940+941+942+819+947+948+822+949+952+953+955+957+834+965+966+847+848+849+850+975+852+977+854+981+984+985+988+861+989+863+991+865+866+867+992+995+871+999+1002+877+1006+1007+880+881+1010+1011+887+889+538+539+540+541+542+543+544+547+550+551+552+553+556+557+559+560+562+565+566+567+569+574+577+586+597+598+599+601+602+603+605+606+607+609+611+623+625+627+628+629+645+648+651+653+654+659+660+661+693+695+698+700+702+703+704+705+706+707+709+710+711+713+714+715+716+717+719+721+722+733+744+747+748+750+759+760+762+775+777+778+781+784+785+787+788+791+792+795+798+799+801+802+803+805+806
PyMol> sele exosc8_con, resid 547+548+551+552+553+559+606+651+653+654+721+760+775+776+777+778+780+781+784
PyMol> color red, sas
PyMol> color magenta, exosc8_con
PyMol> show sticks, exosc8_con
```

<!-- This residues correspond to `16, 17, 20, 22, 75, 120, 122, 123, 190, 229, 230, 244, 245, 246, 247, 249, 250, 253` in the original PDB -->


## Final docking


Run another docking using the filtered residues of the EXOSC8 and Nsp8
* MolA = Nsp8 - Hexadecamer - All surface solvent accessible passive
* MolB = EXOSC8 - 2nn6 (renumbered) - Active: 585,586,589,591,644,689,691,692,759,798,799,813,814,815,816,818,819,822

https://bianca.science.uu.nl/haddock2.4/run/1524990235/31391-nsp8-hexadeca-exosome

Again using "Optimize for Bioinformatics prediction" option

* [Docking parameter](runs/final-docking.json)

```
Hexadecamer Nsp8 solvent accessible:
2304, 2306, 2309, 2310, 2313, 2317, 2328, 2329, 2330, 2332, 2341, 2344, 2345, 2348, 2349, 2350, 2353, 2356, 2358, 2360, 2362, 2363, 2368, 2370, 2371, 2373, 2375, 2376, 2378, 2379, 2380, 2381, 2383, 2384, 2386, 2285, 2286, 2287, 2290, 2294, 2301, 2302

Full exosome target residues
585,586,589,591,644,689,691,692,759,798,799,813,814,815,816,818,819,822
```

## Selection

From the run below we select the Top 10 clusters, alternating between the best model of each cluster and the second, then we fill the selection to 100 models by adding conformations based on their single-structure ranking.

The following script will:

Read the cluster information, rank the clusters and identify to which cluster each structure belongs to. For prosperity, it prepares a dataframe that can easily be parsed with pandas, example: 
```python
>>> df = pd.DataFrame(data, columns=['pdb', 'single_structure_ranking', 'overall_cluster_ranking', 'internal_cluster_ranking'])
>>> top1_top10_cluster = df[(df['overall_cluster_ranking'] <= 5) & (df['internal_cluster_ranking'] == 1)]['pdb']
```


After the models have been selected, the script then does sequence alignment between each of the models and the template to create a reference numbering dictionary so that the final submission matches the provided template. The renumbered models are made into an ensemble and prepared for submission using pdb-tools's `pdb_mkensemble` and `pdb_tidy`

```
$ python scripts/prepare_submission.py -h
usage: prepare_submission.py [-h] run_path template

positional arguments:
  run_path    Location of the run
  template    Template provided by CAPRI

optional arguments:
  -h, --help  show this help message and exit

$ python prepare_submission.py capri_51_183.brk runs/31349-nsp8-hexadeca-exosome
```

*** 
