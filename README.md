# BioPerf: A Benchmark Suite to Evaluate High-Performance Computer Architecture on Bioinformatics Applications

## Motivation:

The exponential growth in the amount of genomic data has spurred growing interest in large scale analysis of genetic information. Bioinformatics applications, which explore computational methods to allow researchers to sift through the massive biological data and extract useful information, are becoming increasingly important computer workloads. This paper presents BioPerf, a benchmark suite of representative bioinformatics applications to facilitate the design and evaluation of high-performance computer architectures for these emerging workloads. Currently, the BioPerf suite contains codes from 10 highly popular bioinformatics packages and covers the major fields of study in computational biology such as sequence comparison, phylogenetic reconstruction, protein structure prediction, and sequence homology & gene finding. We demonstrate the use of BioPerf by providing simulation points of pre-compiled Alpha binaries and with a performance study on IBM Power using IBM Mambo simulations cross-compared with Apple G5 executions.  

The BioPerf suite (available from [**www.bioperf.org**](http://www.bioperf.org/)) includes benchmark source code, input datasets of various sizes, and information for compiling and using the benchmarks. Our benchmark suite includes parallel codes where available.

## Benchmark Developers:

*   [David A. Bader](http://davidbader.net), NJIT
*   Yue Li, Univ. of Florida
*   [Tao Li](http://www.taoli.ece.ufl.edu/), Univ. of Florida
*   [Vipin Sachdeva](https://www.linkedin.com/in/vsachdeva/)

## BioPerf Papers and Presentations:

*   D.A. Bader, Y. Li, T. Li, V. Sachdeva, ``[BioPerf: A Benchmark Suite to Evaluate High-Performance Computer Architecture on Bioinformatics Applications](https://davidbader.net/publication/2005-byltls/),'' [ _The IEEE International Symposium on Workload Characterization (IISWC 2005)_](https://web.archive.org/web/20190506135227/http://www.iiswc.org/), Austin, TX, October 6-8, 2005\. [Powerpoint Slides](http://github.org/Bader-Research/BioPerf/BioPerf-IISWC2005)
*   D.A. Bader and V. Sachdeva, [Incorporating Life Science Applications into the Architectural Optimizations of Next-Generation Petaflops Systems](BioPerf-gatech-051028.ppt), (poster), Computational Science & Engineering Division Open House, College of Computing, Georgia Institute of Technology, 28 October 2005.
*   D.A. Bader and V. Sachdeva, [Incorporating life sciences applications in the architectural optimizations of next-generation petaflop-system](https://web.archive.org/web/20190506135227/http://www.cc.gatech.edu/~bader/papers/CSB2005.html), The 4th IEEE Computational Systems Bioinformatics Conference (CSB 2005), Stanford University, CA, August 8-11, 2005.
*   D.A. Bader and V. Sachdeva, [BioSPLASH: A sample workload from bioinformatics and computational biology for optimizing next-generation high-performance computer systems](https://web.archive.org/web/20190506135227/http://www.cc.gatech.edu/~bader/papers/BioSplash.html), (Poster Session), 13th Annual International Conference on Intelligent Systems for Molecular Biology (ISMB 2005), Detroit, MI, June 25-29, 2005.
*   D.A. Bader, `` [An Open Benchmark Suite for Evaluating Computer Architecture on Bioinformatics and Life Science Applications](https://web.archive.org/web/20190506135227/http://www.bioperf.org:80/Bader-SPEC2006.pdf),'' [_Proc. SPEC Benchmark Workshop 2006_](https://web.archive.org/web/20190506135227/http://www.spec.org/workshops/2006/), Austin, TX, January 2006. [Powerpoint Slides](http://github.org/Bader-Research/BioPerf/SPEC2006-BioPerf.ppt)

## [Over 100 Studies Cite BioPerf](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=6998941502998479924)

## How to Get the BioPerf Benchmark Suite:

1.  Download the BioPerf benchmark suite from this GitHub repository
2.  Download the Inputs (LARGE FILES) from [Dropbox](https://www.dropbox.com/s/mhfh1f9h75eeuiv/BioPerf-Inputs.tar.gz?dl=0)

### Packages included in BioPerf:

*   **BLAST** (blastp, blastn)
*   **CLUSTALW** (clustalw, clustalw_smp)
*   **FASTA** (fasta34_t, ssearch34_t)
*   **GLIMMER** (glimmer, glimmer-package)
*   **GRAPPA** (grappa)
*   **HMMER** (hmmsearch, hmmpfam)
*   **PHYLIP** (dnapenny, promlk)
*   **TCOFFEE** (tcoffee)

## **BioPerf Input Datasets**

<table border="1" <caption="">

<tbody>

<tr>

<th>Package Name</th>

<th>Executables</th>

<th>class-A</th>

<th>class-B</th>

<th>class-C</th>

</tr>

<tr>

<th rowspan="2">BLAST</th>

<th rowspan="2">blastn  

blastp</th>

<td>Sequence of homo sapiens hereditary haemochromatosis</td>

<th>---</th>

<th>---</th>

</tr>

<tr>

<td>E.Coli sequence against Drosoph database</td>

<td>5 sequence of average length 7500</td>

<td>Input dataset of 20 sequences of about 7000 residues each against the Swissprot database</td>

</tr>

<tr>

<th>CE</th>

<th>ce</th>

<td>1hba.pdb and 4hhb.pdb are different  
types of hemoglobin which is used to transport oxygen</td>

<th>---</th>

<th>---</th>

</tr>

<tr>

<th>CLUSTALW</th>

<th>clustalw</th>

<td>**1a02J_1hjbA:** 50 sequences of length almost 60</td>

<td>**1290.seq:** 66 sequences of length almost 1100</td>

<td>**6000.seq:**320 sequences of length almost 1000</td>

</tr>

<tr>

<th rowspan="2">FASTA</th>

<th rowspan="2">fasta34_t  

ssearch34_t</th>

<td>qrhuld.aa is a query file that contains the human LDL receptor precursor</td>

<th>---</th>

<th>---</th>

</tr>

<tr>

<td>1abseq.pep of 40 residues aligned against 5 sequences of length almost 60</td>

<td>Bacteria genomes of length almost 360,000 base pairs</td>

<td>Bacteria genome of 700,000 aligned with 940,000 base pairs</td>

</tr>

<tr>

<th rowspan="2">GLIMMER</th>

<th rowspan="2">glimmer  

glimmer-package</th>

<td>Haemophilus_influenzae genome of about 4.7 million with glimmer.icm which is a collection of Markov models</td>

<th>---</th>

<th>---</th>

</tr>

<tr>

<td>NC_004061.fna of about 650,000 residues</td>

<td>Bacteria genomes of length almost 2.8 million base pairs</td>

<td>Bacteria genome of about 9 million base pairs</td>

</tr>

<tr>

<th>GRAPPA</th>

<th>grappa</th>

<td>10 bluebell flower species</td>

<td>12 bluebell flower species</td>

<td>13 bluebell flower species</td>

</tr>

<tr>

<th rowspan="2">HMMER</th>

<th rowspan="2">hmmsearch  

hmmpfam</th>

<td>Brine shrimp globin compared against 50 HMM's</td>

<th>---</th>

<th>---</th>

</tr>

<tr>

<td>Sequence of 450 base pairs against database of 500 HMM's</td>

<td>Sequence of 9000 residues compared against 2000 HMM's</td>

<td>Sequence of 9000 residues compared against Pfam database (600M - 7500 HMM's)</td>

</tr>

<tr>

<th rowspan="2">PHYLIP</th>

<th rowspan="2">dnapenny  

promlk</th>

<td>16 Aligned sequences of 2581 characters</td>

<td>---</td>

<td>---</td>

</tr>

<tr>

<td>6 aligned sequences of 39 characters</td>

<td>14 aligned sequences of 230 characters each</td>

<td>92 aligned sequences of almost 220 characters</td>

</tr>

<tr>

<th>PREDATOR</th>

<th>predator</th>

<td>Single aligned sequence of 9100 residues</td>

<td>5 sequences each of length 7500</td>

<td>19 sequences each of length about 7000 residues</td>

</tr>

<tr>

<th>TCOFFEE</th>

<th>tcoffee</th>

<td>5 sequences each of length 250</td>

<td>50 sequences each of length almost 280</td>

<td>50 sequences of average length almost 850</td>

</tr>

</tbody>

</table>

## Installation of BioPerf

*   Unzip <big><tt>BioPerf.zip</tt></big>. <big><tt>BioPerf.zip</tt></big> on unzipping contains 2 files: <big><tt>BioPerf-codes.tar.gz</tt></big> and <big><tt>install-BioPerf.sh</tt></big>. <big><tt>install-BioPerf.sh</tt></big> is the script which installs BioPerf on the user machine.
*   Run <big><tt>install-BioPerf.sh</tt></big>. The user is prompted for the directory for the installation of BioPerf, the default directory of BioPerf is <big><tt>$HOME/BioPerf</tt></big>, else the user input is taken.
*   The script then gunzips and untars the <big><tt>BioPerf-codes.tar.gz</tt></big> into the BioPerf installation directory: <big><tt>BioPerf-codes.tar.gz</tt></big> is the complete package of BioPerf consisting of source codes, pre-compiled executables, installation, compiling and running scripts, and the input datasets of varying sizes for each of the codes.
*   On successful installation, the installation directory of BioPerf is written to <big><tt>$HOME/.bioperf</tt></big> which then exports the environment variable <big><tt>$BIOPERF</tt></big>. Every script used for running executables of BioPerf sources this file <big><tt>$HOME/.bioperf</tt></big>, goes to the appropriate subdirectory of <big><tt>$BIOPERF</tt></big>, and then runs the executables. If this file is moved or the installation directory is moved or renamed, the script will fail to run asking the user to edit the <big><tt>.bioperf</tt></big> file.

## BioPerf directory structure

The installation directory of BioPerf contain the following directories:

*   <big><tt>Binaries</tt></big> : Directory containing pre-compiled x86, PowerPC and Alpha binaries, these binaries are contained in separate sub-directories, Alpha-binaries, x86-binaries (Linux) and PowerPC-binaries (Mac OS). x86 and the PowerPC binaries are included for all the executables, while Alpha binaries are not fully included for all the codes. The subdirectories for each of the platform further have directories for each of the packages.  

*   <big><tt>Inputs</tt></big> : Directory containing all the inputs for each of the executables. There are sub-directories for each of the packages, the inputs for each executable are further categorized into class-A, class-B and class-C based on the sizes of the inputs. Some of the input directories have only one class of input, in which case there are no further subdirectories. This directory is large and available via a separate [Dropbox download](https://www.dropbox.com/s/mhfh1f9h75eeuiv/BioPerf-Inputs.tar.gz?dl=0). In case an attempt is made to run a script for a executable which uses any of these databases, the scripts will look for the databases on the host machine in a directory represented by the environment variable <big><tt>$DATABASES</tt></big>. The databases can be downloaded from this website; if the databases are not found, the script will fail. If the databases are downloaded, then the user needs to set an environment variable <big><tt>$DATABASES</tt></big> to the directory where these databases have been downloaded. The run script will then be able to run.  

*   <big><tt>Outputs</tt></big> : All the scripts in BioPerf are set up to store their output in this directory. The Outputs directory itself has sub-directories by the name of each of the packages.  

*   <big><tt>Scripts</tt></big> : This directory as it's name implies, contains all the scripts used in BioPerf except the use-bioperf.sh which is the wrapper for all these scripts. It has further subdirectories that includes the scripts for running each of the codes for small, medium and large datasets, scripts for compiling each of the codes , running the BioPerf suite and installing BioPerf on your architecture incase your architecture's binaries are not part of the BioPerf package. The script has separate subdirectories for each of the tasks, but the naming itself is fairly intuitive.  

*   <big><tt>Simpoints</tt></big> : This directory contains the simulation points using the Simpoint methodology to simulate the representative workload execution phases.  

*   <big><tt>Source-codes</tt></big> : This directory contains the source-codes for each of the codes in BioPerf. The directory structure is same with sub-directories for every package, and further with the executable for the packages having more than one executable.  

    The following scripts are included in BioPerf in the <big><tt>Scripts</tt></big> directory:  

    *   <big><tt>use-bioperf.sh</tt></big> : The wrapper for scripts <big><tt>install-codes.sh</tt></big>, <big><tt>display-versions.sh</tt></big>, <big><tt>run-bioperf.sh</tt></big> and <big><tt>CleanOutputs.sh</tt></big>. These scripts are explained below.  

    *   <big><tt>install-BioPerf.sh</tt> </big>: installs the BioPerf into the directory specified by the user (explained in INSTALLATION section).  

    *   <big><tt>install-codes.sh</tt></big> : compiles the codes for your architecture. Incase your architecture is not x86 with Linux OS, PowerPC with Mac OS or Alpha, you can use this script to compile the codes. This script picks up the makefiles of the source codes from the Source-codes directory, tries to do a make for each of the codes and installs the compiled codes into a subdirectory called <big><tt>$HOSTNAME/Binaries</tt></big>. It will also create <big><tt>$HOSTNAME-scripts</tt></big> subdirectory inside the directory <big><tt>Scripts</tt></big>, which can then be used to run the newly compiled executables. Incase the codes cannot be compiled on your architecture, the script will output an error message telling the user which code failed to compile.  

    *   <big><tt>display-versions.sh</tt></big> : This script outputs the versions of all the installed codes in BioPerf.  

    *   <big><tt>run-bioperf.sh</tt></big> : This script is the basic run script for BioPerf. Its working is explained below in a separate section.  

## How to Use BioPerf

On successful installation, run the script <big><tt>use-bioperf.sh</tt></big> located in the main directory of BioPerf suite to run all the supported tasks in BioPerf.  

The following choices are available:  

*   Run BioPerf [R]
*   Install BioPerf on the user architecture [I]
*   Clean Outputs in the user's <big><tt>$BIOPERF/Outputs</tt></big> [C]
*   Display versions of all installed codes [D]  

*   **Run BioPerf**: If the user selects the option [R], the BioPerf suite is run through the script <big><tt>$BIOPERF/Scripts/Run-scripts/run-bioperf.sh</tt></big>. The script prompts the user for choosing either the platform, if the platform is x86, ppc or alpha, it then gives two modes of running BioPerf: either the user can run all the codes one-by-one, or the user can choose codes which would be run. The user can then add packages to be run, with a prompt for every package. Incase a package is selected to be run, all the executables inside the package are run. After the packages have been selected for running (either all of them if the user chooses the first option or adds the packages one-by-one), the user is prompted for the size of the input datasets which it would like to run for. The size of the input that is selected, all the codes are run for the same input dataset for e.g. if the user selects class-C, all the packages selected for running are run for class-C input dataset. As explained above, some of the larger databases are not included in the BioPerf package to reduce the size of the package as much as possible. If the user selects packages with input sizes (class-C) which require these databases, and the databases are not installed or the <big><tt>$DATABASES</tt></big> environment variable is not set to the appropriate directory, the user is notified, and is then given the choice to run BioPerf without the particular package, or abort the whole running. When BioPerf is subsequently run, it outputs time of running for each executable and also the complete execution time for running the suite with the selected packages.  

*   **Install BioPerf on the user architecture**: If the user selects the option [I], user architecture's binaries are compiled and installed through the script <big><tt>$BIOPERF/Scripts/Run-scripts/install-bioperf.sh</tt></big>. All the codes are attempted to be compiled, and the script fails incase any of the codes fails to compile. The script also tries to first delete the previous installation if detected before trying to proceed with the new installation. This is done because the scripts' first step is to make 2 directories, <big><tt>$HOSTNAME-Binaries</tt></big> and <big><tt>$HOSTNAME-scripts</tt></big> inside the <big><tt>Binaries</tt></big> and the <big><tt>Scripts</tt></big> directory of the BioPerf package respectively. In case of successful compilations for each of the packages, the same sub-directory structure is maintained as in the x86 and the PowerPc sub-directories also. This allows the main <big><tt>run-bioperf.sh</tt></big> to use these executables just as they use x86 and the PowerPC executables.  

*   **Clean Outputs**: If the user selects the option [C], the outputs in <big><tt>$BIOPERF/Outputs</tt></big> are deleted through the script <big><tt>$BIOPERF/Scripts/Run-Scripts/CleanOutputs.sh</tt></big>. All the scripts in BioPerf store the outputs generated in running the executables in the <big><tt>Outputs</tt>></big> directory, which has sub-directories for each of the packages also. This increases the size of the <big><tt>Outputs</tt></big> directory, and hence the script deletes all the outputs previously generated.  

*   **Display Versions**: If the user selects the option [D], the versions of all the software in the BioPerf suite are displayed through the script <big><tt>$BIOPERF/Scripts/Run-scripts/display-versions.sh</tt></big>.


