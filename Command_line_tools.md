## Here are some useful command line tools for working with sequence data at the command line. 

If you want to save output to a new file, use a ">":

```
head copy2.txt > copy2head.txt
```

If you want to search for a term within a txt file, use ```grep```:

```
grep Innexin copy2.txt
```

This will print the output to the terminal. If you want to save the output to a file, use the ">":

```
grep innexin copy2.txt >copy2innexin.txt
```

This is a sequence file in FASTA format. These ```grep``` searches are only returning the file header. What if you want the sequence? Add "-A" for all lines that follow, and for FASTA, add "1" since the following line contains the sequence associated with that header.

```
grep -A 1 innexin copy2.txt
```

```grep``` searches are sensitive to capitalization, etc. To make it less sensitive, use -i

```
grep -A 1 -i Innexin copy2.txt
```

If you want to copy that output to a new file:

```
grep -A 1 -i Innexin copy2.txt > copy2output.txt
```

If you want counts of the nucleotides/amino acids in a sequence, use the ```wc``` word count command:

```
grep -A 1 -i Innexin copy2.txt | wc
```

or gives you the number of hits:

```
grep -A 1 -i -c Innexin copy2.txt
```
If you want to see a history of your commands, and if you want to save the output of that history, run these:

```
history
history > history.txt
```

If you need help, try [Explain Shell](http://explainshell.com/)
