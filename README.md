cgMLSTFinder
===================

Core genome Multi-Locus Sequence Typing
cgMLSTFinder runs KMA [1] against a chosen core genome MLST (cgMLST) database and outputs the detected alleles in a matrix file.


Documentation
=============

The cgMLSTFinder service contains one python script cgMLST.py which is the script of the latest version of 
the cgMLSTFinder service. The method assigns core genome MLST alleles based on WGS data for several isolates. If 
3 or more isolates are run together in one go, and the '--neighbor' flag is set with a executable neighbor 
program, the allele dictances between the isolates will be calculated and a neighbor joining tree will be generated.

## Content of the repository
1. cgMLST.py     - the program
2. README.md
3. Dockerfile   - dockerfile for building the cgmlst docker container
4. make_nj_tree.py - neighbor program


## Installation

Setting up cgMLSTFinder program
```bash
# Go to wanted location for cgmlstfinder
cd /path/to/some/dir
# Clone and enter the cgmlstfinder directory
git clone https://bitbucket.org/genomicepidemiology/cgmlstfinder.git

cd cgmlstfinder
```

Build Docker container
```bash
# Build container
docker build -t cgmlstfinder .
```

Download and install cgMLST database
```bash
# Go to the directory where you want to store the cgmlst database
cd /path/to/some/dir
# Clone install script from git repository
git clone https://bitbucket.org/genomicepidemiology/cgmlstfinder_db.git
cd cgmlstfinder_db
cgMLST_DB=$(pwd)
# Install cgMLST database (look at https://bitbucket.org/genomicepidemiology/cgmlstfinder_db.git for more information)
python3 INSTALL.py
```

This script will install the already kma indexed cgMLST schemes. 

## Dependencies
In order to run the program without using docker, Python 3.5 (or newer) should be installed along with the following versions of the modules (or newer).

#### Python module
- ete3

Module can be installed using the following command. Here, the installation of the module ete3 is used as an example:
```bash
pip3 install ete3
```
#### KMA
KMA should be installed. The instructions here will install the newest version of KMA in the default location cgMLST uses. 
KMA can be installed in another location but the path to KMA will then need to be specified every time you run cgMLSTFinder unless you add the kma program to your PATH.
```bash
# Go to the directoy in which you installed the cgMLSTFinder tool
cd /path/to/some/dir/cgmlstfinder
git clone https://bitbucket.org/genomicepidemiology/kma.git
cd kma && make
```
Further information can be found here:
```url
https://bitbucket.org/genomicepidemiology/kma
```

#### Prodigal (optional)
Make sure prodigal is installed when running preassembled partial, complete genomes.
```url
https://github.com/hyattpd/Prodigal/wiki/installation
```


## Usage
Run Docker container

```bash
# Run cgMLST container
docker run --rm -it \
       -v $cgMLST_DB:/database \
       -v $(pwd):/workdir \
       cgmlstfinder -o [OUTPUT PATH] -s [SPECIE] -db [DATABASE PATH] -t [TEMPORARY FILE] [INPUT/S FASTQ]
```

The program can be invoked with the -h option to get help and more information of the service.

```bash
usage: cgMLST.py [-h] -i INPUT_FILE(S) -s SPECIES [-db DB_DIR] [-o OUTPUT_DIR] [-t TMP_DIR]
                 [-k KMA_PATH] [-p PRODIGAL_PATH] [-n NJ_PATH] [-mem] [-q]
                 


  -i INPUT --inputfile        FASTQ file(s) or FASTA file to do cgMLST on. 
						List several filenames in a comma separeted manner without white spaces 
  -s SPECIES, --species SPECIES
                        Species. Must match the name of a species in the
                        database
optional arguments:
  -h, --help            show this help message and exit

  -db DB_DIR, --databases DB_DIR
                        Directory containing the databases and gene lists for
                        each species.
  -o OUTPUT_DIR, --outdir OUTPUT_DIR
                        Output directory.
  -t TMP_DIR, --tmp_dir TMP_DIR
                        Temporary directory for storage of the results from
                        the external software.
  -k KMA_PATH, --kmapath KMA_PATH
                        Path to executable kma program.
  -p PRODIGAL_PATH --prodigal_path
						Path to Prodigal program
  -n NJ_PATH, --nj_path NJ_PATH
                        Path to executable neighbor joining program.
  -mem, --shared_memory
                        Use shared memory to load database.
  -q --quiet

```

Examples of command to run cgMLSTFinder:

```bash
python3 cgMLST.py /path/to/isolate.fq.gz -s salmonella -o /path/to/outdir -db /path/to/cgmlstfinder_db/
-k /usr/local/bin/kma -n /usr/local/bin/neighbor

python 3cgMLST.py -i /path/to/isolate_1.fa.gz,/path/to/isolate_2.fa.gz -s salmonella -db /path/to/cgmlstfinder_db/ -o /path/to/outdir
```

## Web-server

A webserver implementing the methods is available at the [CGE website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/cgMLSTFinder/

Citation
=======

Please cite the relevant typing scheme, if you intend to publish your results:

**Campylobacter jejuni/coli (pubMLST)**: Campylobacter Multi Locus Sequence Typing website sited at the University of Oxford (Jolley & Maiden 2010, BMC Bioinformatics, 11:595. publication).

**Clostridium (EnteroBase)**:Global gen, et al.,Global genomic population structure of Clostridioides difficile (2019), bioRxiv Aug. 6, 2019.Publication.

**Escherichia coli (EnteroBase)**: Zhou Z., Alikhan NF, Mohamed K, the Agama Study Group, Achtman M (2020), The EnteroBase user's guide, with case studies on Salmonella transmissions, Yersinia pestis phylogeny and Escherichia core genomic diversity. Genome Res. 2020. 30: 138-152. Publication.

**Salmonella (EnteroBase)**:Alikhan NF, Zhou Z,Sergeant MJ, Achtman M. A genomic overview of the population structure of Salmonella. PLoS genetics 2018, 14 (4): e1007261. Publication.

**Yersinia (EnteroBase)**: Zhou Z., Alikhan NF, Mohamed K, the Agama Study Group, Achtman M (2020), The EnteroBase user's guide, with case studies on Salmonella transmissions, Yersinia pestis phylogeny and Escherichia core genomic diversity. Genome Res. 2020. 30: 138-152. Publication.

References
=======

1. Clausen PTLC, Aarestrup FM, Lund O. Rapid and precise alignment of raw reads against redundant databases with KMA. BMC Bioinformatics 2018; 19:307. 

License
=======

Copyright (c) 2014, Ole Lund, Technical University of Denmark
All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=======
```

Get neighbor from:
http://evolution.genetics.washington.edu/phylip/
