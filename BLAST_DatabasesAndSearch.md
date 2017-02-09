##Instructions for building BLAST databases on your local machine using BLAST+ on MacOSX.

BLAST+ is a program that allows you to build custom BLASTable databases from FASTA-formatted sequence assemblies. These assemblies can include genomic or transcriptomic sequence data. For our purposes, we build BLAST+ databases from transcriptome assemblies to probe for genes of interest.

Before we begin, go to this [link to download BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

If the previous link fails, try [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and follow the link for "installers and source code"


Unzip the file to create the BLAST+ directory. Do this manually or on the command line type:
```
cd ~/Downloads
tar -xvzf ncbi-blast-2.6.0+-x64-macosx.tar.gz
```
The directory comes with an annoying name. Rename it something simple like "BLAST+". Do it manually, or on the command line type:
```
mv ncbi-blast-2.6.0+ BLAST+
```

Let's move the directory somewhere more useful, like our Desktop. Do this manually by dragging the folder, or simply type:
```
mv ~/Downloads/BLAST+ ~/Desktop
cd ~/Desktop/BLAST+
```

The BLAST+ program comes with a 'bin' directory that contains all of the necessary executables for running the program (e.g. blastn, blastp, etc.). Now we want to make three new directories for our *own data* to keep ourselves organized:
```
mkdir Input Output Databases
```

We've created three directories with intuitive names. From here, we need only 
...1.) Create a BLASTable database from a FASTA-formatted sequence assembly and 
...2.) Create a list of probes (nucleotide, amino acid, or protein domain sequences) for your gene of interest to be used as a BLAST query.

To build our BLASTable database, place your sequence assembly (in FASTA format) into the "Input" directory. Do this manually, or type:
```
mv path-to-your-file/filename.fasta ~/Desktop/BLAST+/Input
```

Within the main directory there is a subdirectory called "bin". This contains the executable that allows us to build a database and run various BLAST searches:

```
cd bin
```

To make the database, run the `makeblastdb` command: 

```
~/Desktop/Bioinformatics/BLASTplus/bin/makeblastdb 
-in ~/Desktop/Bioinformatics/BLASTplus/Input/SEQUENCE ASSEMBLY FILE 
-dbtype nucl 
-parse_seqids 
-out ~/Desktop/Bioinformatics/BLASTplus/Databases/NAME FOR YOUR DATABASE
```

### Running a BLAST search
Make a txt file with your sequence queries and put this file into your "Input" directory
Run the below commands, but **change paths as necessary** and look at the flags. I have trouble running these commands without putting the explicit paths. 
Make sure to add output file names; here we are exporting data in HTML format! Need to include -html flag.

```
~/Desktop/Bioinformatics/BLASTplus/bin/tblastn 
-query /Users/pvaelli/Desktop/Bioinformatics/BLASTplus/Input/SEQUENCE QUERY FILE
-db /Users/pvaelli/Desktop/Bioinformatics/BLASTplus/Databases/DATABASE FILE
-out /Users/pvaelli/Desktop/Bioinformatics/BLASTplus/Output/OUTPUT FILE.html -html
```

