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

In an anaconda environment, you can run 'conda install numpy scipy pandas biopython; conda install -c conda-forge apsw' to install these Python dependencies.

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
2. Generate any secondary structure prediction files. They should be in PSIPRED format. Possible computational prediction methods include:
    * JPRED (http://www.compbio.dundee.ac.uk/jpred/)
    * Porter (http://distillf.ucd.ie/porterpaleale/)
    * SSPro (http://scratch.proteomics.ics.uci.edu/)
    * PSIPRED (http://bioinf.cs.ucl.ac.uk/psipred/)
    * DeepCNF/RaptorX Property (http://raptorx.uchicago.edu/StructurePropertyPred/predict/)
    * SPIDER3 (http://sparks-lab.org/server/SPIDER3/)
    * DNSS (http://sysbio.rnet.missouri.edu/multicom_toolbox/tools.html)
3. Create fragments with Rosetta using e.g. 'scripts/run_rosetta.py -l <IDP chain ID> -s <IDP fasta file> --psipred_path <psipred ss prediction file> --porter_path <porter ss prediction file>  --jpred_path <jpred ss prediction file> --sspro_path <sspro ss prediction file> -d ./ --nfrag 30 <PDB ID>'
4. Convert Rosetta output to PDB format using e.g. 'scripts/rosetta_to_pdb.py -p <PDB ID> -l <IDP chain ID> -s <IDP fasta file> -f output_files/<PDB ID><IDP chain ID>.30.9mers'
5. Rename 'window_data.csv' to '${PDBID}_data.csv'
6. Convert CA to full-atom backbones using Pulchra like 'find ./ -name "frag_???.pdb" -print0 | xargs -0 -n1 pulchra'
7. Build side-chains using chosen side-chain modeling software
8. Put fragments into subdirectories as in 'test/4ah2'
9. Run LZerD on each fragment (See http://kiharalab.org/proteindocking/lzerd.php). ZDOCK decoys can be used instead, provided the output files 'complex.1', etc. are renamed like 'model1.pdb', for example using the Perl rename command 'rename "s/complex\./model/g; s/$/.pdb/g"'
10. Score the docked models using GOAP and ITScorePro
11. Copy and edit 'test/test_decoys.sh' for new protein and run as for test protein
