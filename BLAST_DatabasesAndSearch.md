##Instructions for building BLAST databases on your local machine using BLAST+ on MacOSX.

BLAST+ is a command line tool that allows you to build custom BLASTable databases from FASTA-formatted sequence assemblies. These assemblies can include genomic or transcriptomic sequence data. For our purposes, we build BLAST+ databases from transcriptome assemblies to probe for genes of interest.

Note: these instructions are meant for MacOSX users, but conceptually it will be the same across platforms. On the Mac, these tools are run using the terminal.

Before we begin, go to this [link to download BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and follow the link for "Installers and source code..". Download the appropriate version for your computer.

***
###Getting started  

Open up the terminal program on your computer.

First, unzip the BLAST+ file to create the BLAST+ directory. You can do this manually, or in the terminal window paste the following code:
```
cd ~/Downloads
tar -xvzf ncbi-blast-2.6.0+-x64-macosx.tar.gz
```  
The directory comes with an annoying name. Rename it something simple like "BLAST+". Do it manually or paste the following into the terminal:
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

**
###Building a BLASTable database  

To build our database, first place your FASTA-formatted sequence assembly in the "Input" directory. Do this manually, or:
```
mv path-to-your-file/filename.fasta ~/Desktop/BLAST+/Input
```

To make the database, we run the `makeblastdb` command. Copy the code below into a text editor such as [TextWrangler]() and modify the paths and filenames to match your own computer. For me, the code looks like this. You need only change the file paths and file names.  
Note:
+Putting /User/pvaelli/Desktop oddly did not work, but using a "~" sign to represent /Users/pvaelli does work. 
+You will need to remove the hard returns before each flag in the code. I added them to make the code easier to read.

Tidy version:
```
~/Desktop/BLAST+/bin/makeblastdb 
-in ~/Desktop/BLAST+/Input/newt.Trinity.fasta 
-dbtype nucl 
-parse_seqids 
-out ~/Desktop/BLAST+/Database/newt_transcriptome
```
Actual version:
```
~/Desktop/BLAST+/bin/makeblastdb -in ~/Desktop/BLAST+/Input/newt.Trinity.fasta -dbtype nucl -parse_seqids -out ~/Desktop/BLAST+/Database/newt_transcriptome
```
What is this code doing? In the first line, we give the path to the executable file 'makeblastdb' in the bin directory. '-in' is a flag that represents input file. '-dbtype' is a flag that tells the program if your FASTA data is nucleotide or amino acid data. '-out' tells the program where to save your new database and what to call it. 

**
### Running a BLAST search
Make a txt file with your sequence queries and put this file into your "Input" directory
Run the below commands, but **change paths as necessary** and look at the flags. I have trouble running these commands without putting the explicit paths. 
Make sure to add output file names; here we are exporting data in HTML format! Need to include -html flag.

Tidy version:
```
~/Desktop/BLAST+/bin/tblastn 
-query ~/Desktop/BLAST+/Input/Homo_nav1.6.fasta 
-db ~/Desktop/BLAST+/Database/newt_transcriptome 
-out ~/Desktop/BLAST+/Output/Nav1.6_hits.html -html
```
Actual version:
```
~/Desktop/BLAST+/bin/tblastn -query ~/Desktop/BLAST+/Input/Homo_nav1.6.fasta -db ~/Desktop/BLAST+/Database/newt_transcriptome -out ~/Desktop/BLAST+/Output/Nav1.6_hits.html -html
```
**
### Specific example

Making a database from transcriptomic data.

Tidy version:
```
~/Desktop/BLAST+/bin/makeblastdb 
-in ~/Desktop/BLAST+/Input/newt.Trinity.fasta 
-dbtype nucl 
-parse_seqids 
-out ~/Desktop/BLAST+/Database/newt_transcriptome
```
Actual version:
```
~/Desktop/BLAST+/bin/makeblastdb -in ~/Desktop/BLAST+/Input/newt.Trinity.fasta -dbtype nucl -parse_seqids -out ~/Desktop/BLAST+/Database/DATABASE_FILE
```

Performing a tBLASTn search for the beta actin gene:

Tidy version:
```
~/Desktop/BLAST+/bin/tblastn 
-query ~/Desktop/BLAST+/Input/Homo_ATCB.fasta 
-db ~/Desktop/BLAST+/Database/newt_transcriptome 
-out ~/Desktop/BLAST+/Output/ATCB_hits.html -html
```
Actual version:
```
~/Desktop/BLAST+/bin/tblastn -query ~/Desktop/BLAST+/Input/Homo_ATCB.fasta -db ~/Desktop/BLAST+/Database/newt_transcriptome -out ~/Desktop/BLAST+/Output/ATCB_hits.html -html
```
