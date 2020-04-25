# Vaccinology Tools

Determine cross-reactivity between peptide sequences based on short sliding windows and/or helical wheel homology of coiled-coil proteins

## Dependencies:

1. Linux or macOS based operating system
2. [EMBOSS program](http://emboss.open-bio.org/html/adm/ch01s01.html)
 - **_NOTE:_** The download instructions above include an FTP site. If you are unable to download from websites via FTP, try [this fedora repo](https://src.fedoraproject.org/lookaside/extras/EMBOSS/) instead.
3. Python3
4. The following Python libraries: os, glob, sys, copy, linecache, subprocess, re, collections, statistics, numpy, biopython, shutil
5. MARCOIL   

## Installation Instructions for MacOS

These instructions are assuming you are compiling and using this software in your home directory, and that you have wget and git installed. If you are installing somewhere else, you will need to adjust the commands accordingly. If you do not have git and wget installed, please check out the documentation to install [homebrew](https://brew.sh/) for macos, as it is the simplest way to get them. 

### 1. Clone this git repo

- `git clone https://github.com/michellearanha/Vaccinology_Tools.git`

### 2. Download and Install EMBOSS.

- Download EMBOSS via [ftp](http://emboss.open-bio.org/html/adm/ch01s01.html) or using the [fedora repository](https://src.fedoraproject.org/lookaside/extras/EMBOSS/)
- Untar the tarball: `tar -xvf EMBOSS-6.6.0.tar.gz`
- In the new EMBOSS directory, run `./configure --prefix=/Users/$('whoami')/EMBOSS-6.6.0/`
- Run `make`
- Run `make install`
- Make some coffee, the last two commands may take a bit to install. 
 
### 3. Install Python3 and required packages.
- Install Python3 if not already installed (https://www.python.org/downloads/)
- Install required Python packages with Pip3 install requirements.txt -r

### 4. Compile MARCOIL
- Unizip the MARCOIL zip file. Navigate to the directory and enter the command 'make'. Then 'make clean'.

 
### 5. Edit your ~/.bash_profile

- Need to include this line: `PATH=$PATH:</directory/containing/Vaccinology_Tools/Heptad_homology>` 
	- Assuming this is in your home directory, the easiest way to do this is to run:

	`echo "PATH=$PATH:~/Vaccinology_Tools/Heptad_homology" >> ~/.bash_profile && source ~/.bash_profile`

### 6. Ensure the heptad_id file is executable

- Check the permissions on the file to make sure you are able to execute it.
	- Easiest way to check is to `cd` to the `~/Vaccinology_Tools/Heptad_homology` directory and run `ls -la`
	- If you need to add the executable flag, just run `chmod +x heptad_id`

## Installation for Linux

#### 1. Clone this git repo

- `git clone https://github.com/michellearanha/Vaccinology_Tools.git`

#### 2. Download and Install EMBOSS

- Download EMBOSS via [ftp](http://emboss.open-bio.org/html/adm/ch01s01.html) or using the [fedora repository](https://src.fedoraproject.org/lookaside/extras/EMBOSS/)
- Install the x11 development libraries: 
	- Debian and Ubuntu:`sudo apt install libx11-dev`
	- Rhel and Centos: `sudo yum install libX11-devel`
	- **_NOTE:_** The next command will complain if you dont install these. 
- Untar the tarball: `tar -xvf EMBOSS-6.6.0.tar.gz`
- In the new EMBOSS directory, run `./configure --prefix=/home/$('whoami')/EMBOSS-6.6.0/`
- Run `make`
- Run `make install`
- Make some coffee, the last two commands may take a bit to install. 

#### 3. Edit your ~/.bashrc

- Need to include this line: `PATH=$PATH:</directory/containing/Vaccinology_Tools/Heptad_homology>` 
	- Assuming this is in your home directory, the easiest way to do this is to run:

		 `echo "PATH=$PATH:~/Vaccinology_Tools/Heptad_homology" >> ~/.bashrc && source ~/.bashrc`

#### 4. Ensure the heptad_id file is executable

- Check the permissions on the file to make sure you are able to execute it.
	- Easiest way to check is to `cd` to the `~/Vaccinology_Tools/Heptad_homology` directory and run `ls -la`
	- If you need to add the executable flag, just run `chmod +x heptad_id`

## Options for heptad_id

`-f [<.fasta>] (seq.fasta) (Input)`

 - File with sequences listed in fasta format

`-r [<.txt>/<.dat>/...] (register.txt) (Input)`

 - File with a list of heptad registers of all sequences (predicted or actual heptad site of the first residue of each sequence

`-c [<type_enum>] (Identity) (Input)`
  
 - Criteria : Identity, Similarity

`-t [<float_percentage>] (45) (Optional)`
  
 - threshold percentage:  sequences that share "criteria" greater than threshold will be saved in a separate folder.

`-E [<dir_path>] (Input)`

 - Directory of EMBOSS binaries

`-A [<type_enum>] (water) (Input)`
  
 - Global alignment - Needleman-Wunsch algorithm : needle

 - Local alignment - Waterman algorithm : water


## Examples 

#### Running a sliding window

`sliding_debug1 [-f <file:fasta file with sequences>] [-w <window size>] [-g <gap size>]  [-t <threshold: percentage>] [-c <String:Criteria-Identity/Similarity>]  [-E <Directory_of_EMBOSS_program>] [-A <alignment_algorithm>] [-l <length of the longest sequence>] [-D <output folder name>] `

#### Running a heptad homology program

`heptad_id [-f <file:fasta file with sequences>][-r <file:heptad register file>][-c <String:Criteria-Identity/Similarity>][-t<threshold: percentage>][-E <Directory_of_EMBOSS_program>] [-A <alignment_algorithm>]`

#### Specific examples:

`sliding_debug1 -f ~/Vaccinology_Tools/Heptad_homology/seq.fasta -w 10 -g 4 -t 25 -c Identity -E ~/EMBOSS-6.6.0/bin/bin -A needle -l 40 -D ~/Vaccinology_Tools/Heptad_homology/NTC6`

`sliding_debug1 -f ~/Vaccinology_Tools/Heptad_homology/seq.fasta -w 15 -g 4 -t 50 -c Similarity -E ~/EMBOSS-6.6.0/bin/bin -A needle -l 40 -D ~/Vaccinology_Tools/Heptad_homology/NTC6`

`heptad_id -f ~/Vaccinology_Tools/Heptad_homology/seq.fasta -r ~/Vaccinology_Tools/Heptad_homology/register.txt -c Identity -t 25 -E ~/EMBOSS-6.6.0/bin -A needle`

`heptad_id -f ~/Vaccinology_Tools/Heptad_homology/seq.fasta -r ~/Vaccinology_Tools/Heptad_homology/register.txt -c Similarity -t 60 -E ~/EMBOSS-6.6.0/bin -A water`


## Definitions and Long Explanations

What the program does:

- Sliding window:

Scanning for matches between sliding k-mer sequence of vaccine type and entire sequence of non-vaccine types (Sliding window approach).

Consider a set of potentially cross reactive protein sequences that are related within species (e.g. antigenic M proteins of different strains of S.pyogenes ) or between species (e.g. antigenic proteins of different viruses belonging to the Flaviviridae virus family such as Dengue, Zika, West Nile that are known induce cross reactive antibodies). Each sequence in the given set is iteratively made the reference sequence and the tool uses a sliding window approach to calculate a pairwise alignment and the number of residues that are identical between successive fragments of the reference sequence and the sequence of the remaining proteins in the set. Matches in terms of sequence similarity as a criterion can also be chosen instead of sequence identity. The subsequence size/ window length/ fragment size (k) for the reference sequence can be adjusted, and the number of amino acids by which the window moves down the reference sequence can be adjusted by adjusting the gap length (m). This gap also determines the degree of overlap between successive refence sequence fragments. The sliding window direction is from the N-terminal to the C-terminal. The number of identical residues between vaccine and non-vaccine type is plotted as a function of window number. If desired, an overall sequence identity or similarity cutoff can be chosen to limit the number of sequences compared. 

- Helical wheel homology:

A standard or canonical coiled coil structure is built by two or more helices twisting around each other forming bundles with their side chains interlocking in a ‘knobs’ into ‘holes’ packing. The regular meshing of knobs into holes packing requires recurrent positions of the side-chains every seven residues along the helix interface. This seven-residue sequence repeat is called a heptad repeat and the positions in the heptad repeat are labelled a-g. The core forming positions (a and d) are usually occupied by hydrophobic residues whereas the remaining, solvent exposed positions (b, c, e, f and g) are dominated by hydrophilic residues. For proteins with coiled-coil secondary structure, machine learned models implemented in webservers exist that can assign the residues to the heptad pattern seen in coiled-coil proteins. We compared the assignment of individual residues in a coiled coil sequence to the heptad by the program MARCOIL (website) and found that it was identical to the actual heptad assignment of three M protein coiled coil crystal structures (pdb ids: 2oto, 5hyt, 5hzp) that were of interest to us. Similar assignments can also be done using NCOILS, PCOILS. 

The tool requires three inputs 1) a text file with the sequences of all vaccine and non-vaccine types in the fasta format, 2) a text file with identifiers/headers of the fasta sequences without ‘>’ character 3) a file with the heptad registers of all sequences (assigned heptad position of the first residue of all sequences) which can be obtained using MARCOIL. The identity at corresponding heptad positions of each peptide/protein with every other peptide/protein in the given set is output in a tabular format. Additionally, an output is also provided as a matrix in tabular format that contains heptad identity of every sequence with every other sequence which can be clustered to obtain clusters of cross-reactive peptides. An R-based script for clustering is also provided. 


## A step-by-step guide to using all of the tools with an example

1. Navigate to the MARCOIL directory and enter the following command: ./marcoil +dlSs /PATH/TO/FASTA/SEQS/example_seqs.txt
2. Move the three output files (ProbPerState, ProbList, and Domains) located in the Outputs directory to a directory of your choosing for safe-keeping.
3. Navigate to the directory MARCOIL_ANALYSIS. Enter the following command: 

   python3 PATH/TO/PROBLIST/FILE/ProbList ~/Vaccinology_Tools/Heptad_homology/MARCOIL_RESULTS_FOR_INPUT_ANTIGENS/ 35

4. Navigate to your output, specifically, ~/Vaccinology_Tools/Heptad_homology/MARCOIL_RESULTS_FOR_INPUT_ANTIGENS/EQUAL_SIZED_FRAGS_MARCOIL_ANALYSIS. Copy the files CONDENSED_PERFECT_CC_SEQS.txt and REGISTER_FOR_CONDENSED_PERFECT_CC_SEQS.txt to an empty directory, ~/Vaccinology_Tools/Heptad_homology/DIRECTORY_FOR_EXAMPLE_RUN/.

5. Assuming EMBOSS and both the heptad homology and heptad scoring scripts are in your bash profile The following command should run: 

   heptad_id -f ~/Vaccinology_Tools/Heptad_homology/DIRECTORY_FOR_EXAMPLE_RUN/CONDENSED_PERFECT_CC_SEQS.txt-r ~/Vaccinology_Tools/Heptad_homology/DIRECTORY_FOR_EXAMPLE_RUN/REGISTER_FOR_CONDENSED_PERFECT_CC_SEQS.txt -c Identity -t 80 -E ~/EMBOSS-6.6.0/bin -A needle
   
6. Assuming the Ideal coiled coil builder is also in your path enter the following command:
	
   python3 Ideal_CC_Builder.py ~/Vaccinology_Tools/Heptad_homology/DIRECTORY_FOR_EXAMPLE_RUN/CONDENSED_PERFECT_CC_SEQS.txt ~/Vaccinology_Tools/Heptad_homology/DIRECTORY_FOR_EXAMPLE_RUN/REGISTER_FOR_CONDENSED_PERFECT_CC_SEQS.txt 4 10 ~/Vaccinology_Tools/Heptad_homology/DIRECTORY_FOR_EXAMPLE_RUN/IDEAL_COILED_COIL_PROTEINS/
   
## An explanation on the ideal coiled coil building process
 The workflow begins with the user running MARCOIL on an input FASTA file containing the sequences of all antigenic variants that need to be "covered" by the vaccine in question. After the MARCOIL analysis the generated output must be parsed to identify "regular fragments" which are of a user-defined length. Here, regular fragment means MARCOIL has given the fragment a register with no "breaks" in its heptad assignment. For example, a fragment with the assignment "abcdefgabcdefgabcdefg" would be a regular fragment with three heptads. Another example of a regular fragment would be "cedfgabcedfgabcedfgabcedfgab", which has 4 regular heptads. An irregular sequence would not have the heptad periodicity. For example, "cedGdab". There's a break in the periodicity starting at G. For the purposes of vaccine sequences, it is important to remember that the coiled coil fragments represent epitopes. For the purposes of yielding functional antibodies against a particular antigen it is advisable to use 4-5 ideal heptads (28-35 AA). Three heptads (21 AA) can be used if necessary, however more heptads mrean more overlapping epitopes and potentially more cross-reactivity. 
 
 Once the regular heptads are identified, the regular heptads are cut into fragments of equal length. For example, if the argument from step 3 is 35. All regular fragments at this point are greater than or equal to 35 AA in length. Next all of the regular fragments are cut into overlapping 35mers with a sliding window of 1 AA. For example, assume MARCOIL analysis generates the following regular fragment:
 
>emm_2_1-50_regular_fragment_1
VKKEAKLSEAELHDKIKNLEEEKAELFEKLDKVEEEHKKVEEE

The fragment is longer than 35 AA so it will be cut into 35 AA sequences like so:

>emm_2_1-50_regular_fragment_1_window_1
VKKEAKLSEAELHDKIKNLEEEKAELFEKLDKVEE
>emm_2_1-50_regular_fragment_1_window_2
KKEAKLSEAELHDKIKNLEEEKAELFEKLDKVEEE
>emm_2_1-50_regular_fragment_1_window_3
KEAKLSEAELHDKIKNLEEEKAELFEKLDKVEEEH
>emm_2_1-50_regular_fragment_1_window_4
EAKLSEAELHDKIKNLEEEKAELFEKLDKVEEEHK
>emm_2_1-50_regular_fragment_1_window_5
AKLSEAELHDKIKNLEEEKAELFEKLDKVEEEHKK
>emm_2_1-50_regular_fragment_1_window_6
KLSEAELHDKIKNLEEEKAELFEKLDKVEEEHKKV
>emm_2_1-50_regular_fragment_1_window_7
LSEAELHDKIKNLEEEKAELFEKLDKVEEEHKKVE
>emm_2_1-50_regular_fragment_1_window_8
SEAELHDKIKNLEEEKAELFEKLDKVEEEHKKVEE
>emm_2_1-50_regular_fragment_1_window_9
EAELHDKIKNLEEEKAELFEKLDKVEEEHKKVEEE

If the sequence is already exactly 35 AA in length, it is not cut into overlapping k-mers because there are no other options.

Once these all 35mers from the regular fragments are generated, MARCOIL is run a second time to obtain the heptad register for each fragment. Any fragment with an irregular heptad register is pruned. At this point there are 35mer CC sequences with a regular heptad periodicity. It is likely that there are some redundant 35mers derived from different variants of the antigens in your output. For example:

>emm_2.0_1-50_regular_fragment_1_window_1
VKKEAKLSEAELHDKIKNLEEEKAELFEKLDKVEE
>emm_2.1_1-50_regular_fragment_1_window_1
VKKEAKLSEAELHDKIKNLEEEKAELFEKLDKVEE

There is no need to perform alignments for both Emm2.0 and subtype Emm2.1. These will be compiled to a single sequence. For example:
>emm_2.0_1-50_regular_fragment_1_window_1_IDENTICAL_emm_2.1_1-50_regular_fragment_1_window_1
VKKEAKLSEAELHDKIKNLEEEKAELFEKLDKVEE

This dramatically reduces both the size and time the Needleman-Wunsch alignments will take to run. At the end of step 3 you should have two files to use in the Needleman-Wunsch alignments:
1. CONDENSED_PERFECT_CC_SEQS.txt
2. REGISTER_FOR_CONDENSED_PERFECT_CC_SEQS.txt

The Needleman-Wusnch alignments between all fragments in the input FASTA will run calculating linear sequence identity between all of the fragments at a user-defined threshold. The results of these calculations are written to a folder titled HIGH_PERCENTAGE_IDENTITY. The script needleparser.py Parses out all identity percentages stored within each .needle alignment file. These percentages are written in a TSV format to a directory titled PERCENTAGE_IDENTITY_. Next, the data has a column titled "Antigen_Name_Only" appended. Antigen name only makes a new row for each antigen that is stored within the concatenated FASTA header. For example, the concatenated header is:

>emm_2.0_1-50_regular_fragment_1_window_1_IDENTICAL_emm_2.1_1-50_regular_fragment_1_window_1

There are two pairs of delimiters to get the antigens:

Pair 1: > and _regular
Pair 2: _IDENTICAL and _regular

Pair 1 always gets the antigen directly after the FASTA header, in this case that is: emm_2.0_1-50
Pair 2 Gets any other antigen in the FASTA header, in this case that is: emm_2.1_1-50

The column Antigen Name Only allows each percent identity file to be open and quickly count how many antigens a particular window matches at a given percentage identity threshold. This will be extremely useful when it is time to build the ideal CC proteins. 

Next the script heptad_analysis_.py calculates, percentage identities at each heptad site, (a-g), are calcuated. For the purposes of cross-reactivity it is important to keep in mind that results from current data (Aranha 2020) suggests that is identity at the polar sites (b,c,e,f, and g) which correlate with cross-reactivity. It is presumed the a and d sites have minimal solvent exposure as they face the inside of the coiled coil. In turn, this minimal exposure leads to a and d sites presumably having minimal controbutions to antibody binding, if any. The empirical heptad score (Aranha 2020) is calculated to predict immunological relationships with the script heptad_analysis_.py. As a rule of thumb, two fragments with a heptad score >=10.5 may be immunologically related. Keep in mind the heptad score was created for antigens with between 4-5 heptads. Therefore heptad scores >=10.5 for fragments with fewer heptads are less likely to be predictive of an immunological relationship. After calculating heptad scores, each fragment has its own separate file written to a directory titled, HEPTAD_SITE_WITH_SCORES_ALL_WINDOWS_GREATER_THAN_10.5. Each file is named after a particular window and the contents of the file are all windows with heptad scores >=10.5 when their sequence is compared to the particular window in question.

The process of automated coiled coil building requires: 1. CONDENSED_PERFECT_CC_SEQS.txt 2. REGISTER_FOR_CONDENSED_PERFECT_CC_SEQS.txt 3. The number of ideal coiled coil proteins to build. 4. The maximum number of fragments that can be used in your vaccine. 5. Build based on heptad score (Identity at b,c,e,f,g sites only). Linear identity (Identity at all heptad sites a-g) 6. Minimum percent identity required if linear identity is used.

The coiled coil building process is best thought of in phases. The first phase is identifying the most conserved _N_ windows, based on the user-specified criteria. For example, say the user specifies he or she wants 30 windows over 3 proteins based on linear identity at a threshold of 80%. First, CONDENSED_PERFECT_CC_SEQS_.txt is read to create a list of headers. These headers are used to open the proper files stored under PERCENTAGE_IDENTITY. Each file in PERCENTAGE_IDENTITY is read in as a PANDAS dataframe and all antigens (as determined by the column Antigen_Name_Only) which are not at >=80% identity are dropped to create a temporary dataframe. The number of rows in the temporary dataframe corresponds to the number of antigens matched. The name of the antigen file is appended to a master dictionary as a key with its number of matched antigens as a value. After all antigens have been analyzed, Python3 determines the antigen file(key), with the most matched antigens (highest value). 

In some cases, there will be more than one "best" window or fragment. For example, the most number of antigens matched for a particular iteration may be three. In this case all fragments/windows matching three antigens are opened. Consider the following example:

>emm_2.0_1-50_regular_fragment_1_window_1
Could match : Emm2.0 (100%)
		 	  Emm73.0 (80%)
			  Emm124.0 (85%)
			  TOTAL: (265/300)= 88.3% 
>emm_11.0_1-50_regular_fragment_1_window_1
Could match : Emm11.0 (100%) 
			  Emm85.0 (95%)
              Emm158 (90%)
			  TOTAL: (285/300)= 95.0% 

emm_11.0_1-50_regular_fragment_1_window_1, has an average of 95% identity across the three antigens it matches, so it is the selected window for the current iteration. If more than one "best" window also has the highest average one of the "best" windows with the highest average is selected.

After the window is selected the antigens matched by the previous window are added to a drop list. For example, if emm_11.0_1-50_regular_fragment_1_window_1 was selected, Emm11.0, Emm85.0, Emm158.0 would be added to a drop list because they are now represented in the vaccine. Every iteration after the first runs the same way. The match dictionary is cleared. The percentage identity files are opened as PANDAS dataframes and any row which is not >= the user-defined identity threshold is dropped. Additionally, any row which is for an antigen that has previously been covered is now dropped. This ensures the best window available in each iteration is selected.

It should also be noted that if all query antigens are "covered"/represented by the vaccine the program will stop gathering new windows. For example, assume that the program finds that all query antigens are represented using 25 windows at 80% identity. The program will stop looking for antigens and attempt to equally distribute 25 windows among three proteins. This will lead to two proteins with 8 windows and one protein with 9 windows.

After the most conserved windows are identified the heptad register for each of these conserved windows is obtained. Because all windows are of equal size and have regular periodicity, the initial heptad letter and the final heptad letter in a fragment are known. The initial heptad letter is stored in REGISTER_FOR_CONDENSED_PERFECT_CC_SEQS.txt. The final heptad letter is known because the length of the sequences are identical.

The conserved windows are pieced together to create idealized coiled coil proteins. The goal while building these proteins is to trim as few residues from the coiled coil fragments as possible. Residue trimming occurs from the N-terminus and is necessary for maintaining the periodicity of the heptad in the recombinant protein. For example, assume windows of 5 heptads in length are being pieced together to build a recombinant protein.

Assume the previously attached window had an "A" register or was an "A window". The "A window" will have "G" as its heptad letter at amino acid 35. If another window is to be added, the best window to add would be another "A window" because no residue trimming would be necessary. Assume no "A windows" remain. The next best option would be a "G window" with 1 residue trimmed from the N-terminus, then an "F window" with 2 residues trimmed from the N-terminus.

It's not plausible to predict the minimium number of residues that will need to be trimmed to piece together a set of ideal coiled coil windows. So the Python script builds seven candidate vaccine antigens with each iteration. Each antigen starts with a different heptad periodicity (a-g). The vaccine antigen with the fewest residues trimmed is selected as the "best antigen". If there is more than one "best antigen" MARCOIL is run de novo to identify which antigen has the highest average coiled coil probability throughout.

Once the antigen with the fewest residues trimmed and then highest coiled probability is determined the windows/fragments used to build the idealized coiled coil vaccine antigen are dropped from the set of most conserved windows. The windows that remain after the drop are used as candidates in the next iteration. The process repeats itself ultimately building the user-specified number of proteins.

If the user specifies to use heptad score for the building process is identical, but the files under HEPTAD_SITE_WITH_SCORES_ALL_WINDOWS_GREATER_THAN_10.5 are queried for the identification of conserved windows.          