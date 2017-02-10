##Instructions for building BLAST databases on your local computer using BLAST+ on MacOSX.

BLAST+ is a command line tool that allows you to build custom BLASTable databases from FASTA-formatted sequence assemblies. These assemblies can include genomic or transcriptomic sequence data. For our purposes, we build BLAST+ databases from transcriptome assemblies to probe for genes of interest.

Note: these instructions are meant for MacOSX users, but conceptually it will be the same across platforms. On the Mac, these tools are run using the terminal.

Before we begin, go to this [link to download BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and follow the link for "Installers and source code..". Download the appropriate version for your computer. For my Mac, I downloaded ncbi-blast-2.6.0+-x64-macosx.tar.gz

***
###Getting started  

Open up the terminal program on your computer (command + space, search for terminal).

First, unzip the BLAST+ file to create the BLAST+ directory. You can do this manually, or in the terminal window paste the following code (note: change the .tar.gz file name to match the version you downloaded):
```
cd ~/Downloads
tar -xzvf ncbi-blast-2.6.0+-x64-macosx.tar.gz
```  
The directory comes with an annoying name. Rename it something simple like "BLAST+". Do it manually or paste the following into the terminal:
```
mv ncbi-blast-2.6.0+ BLAST+
```  
Let's move the directory somewhere more useful, like our Desktop. Do this manually by moving the folder in Finder, or:
```
mv ~/Downloads/BLAST+ ~/Desktop
cd ~/Desktop/BLAST+
```

The BLAST+ program comes with a 'bin' directory that contains all of the necessary executables for running the program (e.g. makeblastdb, blastn, blastp, etc.). Now we want to make three new directories for our *own data* to keep ourselves organized:
```
mkdir Input Output Databases
```
We've created three directories with intuitive names. From here, we need only:  
1. Create a BLASTable database from a FASTA-formatted sequence assembly and  
2. Create a list of probes (nucleotide, amino acid, or protein domain sequences) for your gene of interest to be used as a BLAST query.

***
###Building a BLASTable database  

To build our database, first place your FASTA-formatted sequence assembly in the "Input" directory. Do this manually, or:
```
mv path-to-your-file/FILENAME.fasta ~/Desktop/BLAST+/Input
```
You will need to adjust the path-to-your-file to match your computer. The assembly was on my Desktop, so:
```
mv ~/Desktop/assembly.fasta ~/Desktop/BLAST+/Input
```

To make the database, we run the `makeblastdb` command. This command has 5 inputs:  
1. The path to the makeblastdb executable file  
2. The path to your input file, which is your FASTA-formatted assembly  
3. Specify if you have nucleotide or amino acid data. After `-dbtype`, add "nucl" or "prot" depending on your file.  
4. `-parse_seqids` tells the program to parse by "|" characters in the sequence headers. You don't need to know any more than that.  
5. The path for the output, which will be our new BLASTable database. Therefore, we will direct the output to the "Database" directory.  

Copy the code below into a text editor such as TextEdit, TextWrangler, or Sublime Text and modify the filenames. If you've followed my instructions thus far, the paths should be the same.

```
~/Desktop/BLAST+/bin/makeblastdb 
-in ~/Desktop/BLAST+/Input/TrinityAssembly.fasta 
-dbtype nucl 
-parse_seqids 
-out ~/Desktop/BLAST+/Database/transcriptome_database
```

Notes on the above code:
* I've added hard returns so you can read the code -- these returns MUST BE DELETED or the code won't work.
* Occasionally, this command can be buggy. If you encounter an error, try replacing the "~" with the full path name. In the case above, it would instead be /Users/pvaelli/Desktop. It will likely be /Users/YOURUSERNAME/Desktop for you. To determine the full path to a folder or file on the terminal in general, navigate to that location using `cd FolderName` and type `pwd`. Navigate backwards up folders using `cd ..`

Example version:
```
~/Desktop/BLAST+/bin/makeblastdb -in ~/Desktop/BLAST+/Input/newt.Trinity.fasta -dbtype nucl -parse_seqids -out ~/Desktop/BLAST+/Database/newt_transcriptome
```

Again, what is this code doing? In the first line, we give the path to the executable file 'makeblastdb' in the bin directory. '-in' is a flag that represents input file. '-dbtype' is a flag that tells the program if your FASTA data is nucleotide or amino acid data. '-out' tells the program where to save your new database and what to call it.  

***
### Running a BLAST search
Now we are almost ready for BLAST. We have made a database that is compatible with BLAST+, but we need to create a FASTA file with query sequence(s) from our gene of interest and put this file into your "Input" directory.  

Here's an example of a FASTA formatted sequence query for the mouse beta-actin (actb) gene:
```
>NM_007393.5 Mus musculus actin, beta (Actb), mRNA
TATAAAACCCGGCGGCGCAACGCGCAGCCACTGTCGAGTCGCGTCCACCCGCGAGCACAGCTTCTTTGCA
GCTCCTTCGTTGCCGGTCCACACCCGCCACCAGTTCGCCATGGATGACGATATCGCTGCGCTGGTCGTCG
ACAACGGCTCCGGCATGTGCAAAGCCGGCTTCGCGGGCGACGATGCTCCCCGGGCTGTATTCCCCTCCAT
CGTGGGCCGCCCTAGGCACCAGGGTGTGATGGTGGGAATGGGTCAGAAGGACTCCTATGTGGGTGACGAG
GCCCAGAGCAAGAGAGGTATCCTGACCCTGAAGTACCCCATTGAACATGGCATTGTTACCAACTGGGACG
ACATGGAGAAGATCTGGCACCACACCTTCTACAATGAGCTGCGTGTGGCCCCTGAGGAGCACCCTGTGCT
GCTCACCGAGGCCCCCCTGAACCCTAAGGCCAACCGTGAAAAGATGACCCAGATCATGTTTGAGACCTTC
AACACCCCAGCCATGTACGTAGCCATCCAGGCTGTGCTGTCCCTGTATGCCTCTGGTCGTACCACAGGCA
TTGTGATGGACTCCGGAGACGGGGTCACCCACACTGTGCCCATCTACGAGGGCTATGCTCTCCCTCACGC
CATCCTGCGTCTGGACCTGGCTGGCCGGGACCTGACAGACTACCTCATGAAGATCCTGACCGAGCGTGGC
TACAGCTTCACCACCACAGCTGAGAGGGAAATCGTGCGTGACATCAAAGAGAAGCTGTGCTATGTTGCTC
TAGACTTCGAGCAGGAGATGGCCACTGCCGCATCCTCTTCCTCCCTGGAGAAGAGCTATGAGCTGCCTGA
CGGCCAGGTCATCACTATTGGCAACGAGCGGTTCCGATGCCCTGAGGCTCTTTTCCAGCCTTCCTTCTTG
GGTATGGAATCCTGTGGCATCCATGAAACTACATTCAATTCCATCATGAAGTGTGACGTTGACATCCGTA
AAGACCTCTATGCCAACACAGTGCTGTCTGGTGGTACCACCATGTACCCAGGCATTGCTGACAGGATGCA
GAAGGAGATTACTGCTCTGGCTCCTAGCACCATGAAGATCAAGATCATTGCTCCTCCTGAGCGCAAGTAC
TCTGTGTGGATCGGTGGCTCCATCCTGGCCTCACTGTCCACCTTCCAGCAGATGTGGATCAGCAAGCAGG
AGTACGATGAGTCCGGCCCCTCCATCGTGCACCGCAAGTGCTTCTAGGCGGACTGTTACTGAGCTGCGTT
TTACACCCTTTCTTTGACAAAACCTAACTTGCGCAGAAAAAAAAAAAATAAGAGACAACATTGGCATGGC
TTTGTTTTTTTAAATTTTTTTTAAAGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAGTTTTTTTGTTTT
GTTTTGGCGCTTTTGACTCAGGATTTAAAAACTGGAACGGTGAAGGCGACAGCAGTTGGTTGGAGCAAAC
ATCCCCCAAAGTTCTACAAATGTGGCTGAGGACTTTGTACATTGTTTTGTTTTTTTTTTTTTTTGGTTTT
GTCTTTTTTTAATAGTCATTCCAAGTATCCATGAAATAAGTGGTTACAGGAAGTCCCTCACCCTCCCAAA
AGCCACCCCCACTCCTAAGAGGAGGATGGTCGCGTCCATGCCCTGAGTCCACCCCGGGGAAGGTGACAGC
ATTGCTTCTGTGTAAATTATGTACTGCAAAAATTTTTTTAAATCTTCCGCCTTAATACTTCATTTTTGTT
TTTAATTTCTGAATGGCCCAGGTCTGAGGCCTCCCTTTTTTTTGTCCCCCCAACTTGATGTATGAAGGCT
TTGGTCTCCCTGGGAGGGGGTTGAGGTGTTGAGGCAGCCAGGGCTGGCCTGTACACTGACTTGAGACCAA
TAAAAGTGCACACCTTACCTTACACAAACAAAAAAAAAAAAAAAA
```

You can paste this sequence into your text editor, and Save As `mouse_actb.fasta` at the correct path location: `~/Desktop/BLAST+/Input`. Make sure to add the extension .fasta to your text file; don't use .txt  

Now we're ready to go. We have a FASTA-formatted query sequence and a BLAST+ friendly database. Here is the code we use to perform a BLASTn (nucleotide query, nucleotide dataset).

Tidy version:
```
~/Desktop/BLAST+/bin/blastn 
-query ~/Desktop/BLAST+/Input/mouse_actb.fasta 
-db ~/Desktop/BLAST+/Database/newt_transcriptome 
-out ~/Desktop/BLAST+/Output/mouse_actb_hits.html -html
```
Actual version:
```
~/Desktop/BLAST+/bin/blastn -query ~/Desktop/BLAST+/Input/mouse_actb.fasta -db ~/Desktop/BLAST+/Database/newt_transcriptome -out ~/Desktop/BLAST+/Output/mouse_actb_hits.html -html
```
So what's going on up here? First, we direct toward the blastn executable in the bin. `-query` specifies the location and name of our query file; `-db` specifies the location and name of our database; and `-out` specifies the location for the output, which is a list of your top BLAST hits. Here, we are creating a new file `mouse_actb_hits.html` using the `-html` flag to create the HTML file. We can then open this file in our web browser to analyze the results of our BLAST search.

Before we wrap up, there are a few different types of BLAST searches that we can perform. Here is a list that includes search type, query type, and data type. By type, I mean either nucleotide data or amino acid data. Here are examples:

| BLAST search  | Query       | Database  |
| ------------- |-------------| ----------|
| blastn        | nucleotide  | nucleotide|
| blastp        | protein     | protein   |
| tblastn       | protein     | nucleotide| 
| blastx        | nucleotide  | protein   |

That's it! You've now gone from a FASTA-formatted sequence assembly to a BLASTable database to a HTML formatted file with a list of your top hits. Now you can copy the name of your top hits and open your original sequence assembly to search and retrieve your hits. You should probably [BLAST your hits on NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to confirm the identity of each sequence. 

Extra:
Here's a [list of all BLAST+ executables](https://www.ncbi.nlm.nih.gov/books/NBK279675/) and their individual commands/flags.
***
### Specific example

Making a database from our transcriptomic data.

Example `makeblastdb`:
```
~/Desktop/BLAST+/bin/makeblastdb 
-in ~/Desktop/BLAST+/Input/newt.Trinity.fasta 
-dbtype nucl 
-parse_seqids 
-out ~/Desktop/BLAST+/Database/newt_transcriptome
```
Remember to remove hard returns. The code should actually look like this:
```
~/Desktop/BLAST+/bin/makeblastdb -in ~/Desktop/BLAST+/Input/newt.Trinity.fasta -dbtype nucl -parse_seqids -out ~/Desktop/BLAST+/Database/DATABASE_FILE
```

Running a BLAST search.

Example `tblastn` searches:
```
~/Desktop/BLAST+/bin/tblastn 
-query ~/Desktop/BLAST+/Input/Homo_ATCB.fasta 
-db ~/Desktop/BLAST+/Database/newt_transcriptome 
-out ~/Desktop/BLAST+/Output/ATCB_hits.html -html
```
or

```
~/Desktop/BLAST+/bin/tblastn -query ~/Desktop/BLAST+/Input/Homo_Kv.fasta -db ~/Desktop/BLAST+/Database/newt_transcriptome -out ~/Desktop/BLAST+/Output/Homo_Kv_hits.html -html
```
