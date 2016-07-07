##Instructions for building BLAST databases on your local machine using BLAST+.

First, go to this [link to download BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

Move this directory (jargon for 'folder') somewhere useful. Maybe your Desktop -- I keep mine in a directory/folder called "Bioinformatics" on my Desktop.

The directory comes with an annoying name. Rename it something simple like "BLASTplus", and move into the directory:

```
cd BLASTplus
```

Now we want to make three new directories to keep ourselves organized:

```
mkdir Input Output Databases
```

We will put our various files into these directories. As you might guess, we will put our assembled sequences/contigs into **Input**, the results of our BLAST searches into **Output**, and our newly created database in, well **Databases**.

To build our BLASTable database, place your sequence assembly (in FASTA format) into the "Input" directory.

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

