IDP-LZerD
=========

IDP-LZerD models the bound conformation of a disordered PPI, where one
intrinsically disordered protein (IDP) binds to an ordered protein.
The inputs are the sequence of a disordered protein and the structure of an
ordered protein.

Copyright (C) 2016-2017 Lenna X. Peterson, Daisuke Kihara, and Purdue University.

License: GPL v3 for academic use.
(For commercial use, please contact us for different licensing)

Contact: Daisuke Kihara (dkihara@purdue.edu)

Reference: Lenna Peterson, A Roy, C Christoffer, G Terashi, and D Kihara. (2017) Modeling disordered protein interactions from biophysical principles. _PLoS Computational Biology_ 13: e1005485. doi: [10.1371/journal.pcbi.1005485](http://dx.doi.org/10.1371/journal.pcbi.1005485)

Installation
============

To run the test, LZerD and all Python dependencies are required.

Python dependencies
-------------------
- apsw
- numpy
- scipy
- pandas
- Biopython

Binary dependencies
-------------------
- blastpgp and nr database (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- LZerD (http://kiharalab.org/proteindocking/lzerd.php)
- Rosetta (https://www.rosettacommons.org)
- Pulchra (http://cssb.biology.gatech.edu/PULCHRA)
- One side-chain modeling software such as:
    * SCCOMP (http://www.sheba-cancer.org.il/cgi-bin/sccomp/sccomp1.cgi)
    * Scwrl4 (http://dunbrack.fccc.edu/scwrl4/)
    * Oscar-star (http://sysimm.ifrec.osaka-u.ac.jp/OSCAR)
    * RASP (http://jianglab.ibp.ac.cn/lims/rasp/rasp)
- GOAP (http://cssb.biology.gatech.edu/GOAP/index.html)
- ITScorePro (http://zoulab.dalton.missouri.edu/resources.html)
- charmm (https://www.charmm.org)

Getting started
===============

0. Edit PATHS.ini to specify path to LZerD, blastpgp, nr, and Rosetta

Test protein
------------
1. Edit 'test/test_decoys.sh' to specify path to IDP-LZerD
2. Run test (NOTE: creates ~250 GB of files): `cd test && ./test_decoys.sh`
3. Non-refined paths will be output in 'test/4ah2ac3'
4. Follow standard procedures to prepare PDB files for CHARMM
5. Edit CHARMM input file 'idp_relax.inp' to run CHARMM

Running a new protein
---------------------
1. Prepare FASTA sequence of IDP and structure of ordered partner protein
2. Create fragments with Rosetta using rosetta_tools/fragment_tools/make_fragments.pl and the flags '-frag_sizes 9 -n_frags 30'
3. Convert Rosetta output to PDB format using 'scripts/rosetta_to_pdb.py'
4. Rename 'window_data.csv' to '${PDBID}_data.csv'
5. Convert CA to full-atom backbones using Pulchra
6. Build side-chains using chosen side-chain modeling software
7. Put fragments into subdirectories as in 'test/4ah2'
8. Run LZerD on each fragment
9. Score LZerD output using GOAP and ITScorePro
10. Copy and edit 'test/test_decoys.sh' for new protein and run as for test protein
