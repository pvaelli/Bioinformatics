# Trinity assembly protocol on Harvard Canon cluster
### based on https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
### revised 10-2-2020

# Step 1: Quality metrics for raw reads
### This job takes about an hour to complete

# submission script:
```
#!/bin/bash
#SBATCH -p shared       # Partition to submit to
#SBATCH -n 1                   # Number of cores
#SBATCH -t 1-8:00               # Runtime in days-hours:minutes
#SBATCH --mem 2000              # Memory in MB
#SBATCH -J FastQC               # job name
#SBATCH -o FastQC.%A.out        # File to which standard out will be written
#SBATCH -e FastQC.%A.err        # File to which standard err will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=<pvaelli@fas.harvard.edu>  # Email to which notifications will be sent

source new-modules.sh
module purge
module load fastqc/0.11.5-fasrc01

for file in *.fastq.gz; do
        fastqc $file
done
```


# Step 2: Identify errorneous kmers

### Script for Rcorrector. Can take up to 3 days

```
#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p shared
#SBATCH -e rcorrect_%A.e
#SBATCH -o rcorrect_%A.o
#SBATCH -J rcorrect
#SBATCH --mem=24000
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=pvaelli@fas.harvard.edu # Email to which notifications will be sent

module purge
module load Rcorrector/20180919-fasrc01 
perl /n/helmod/apps/centos7/Core/Rcorrector/20180919-fasrc01/bin/run_rcorrector.pl -t 12 -1 $1 -2 $2

```

# submit script using this command:
```
sbatch Rcorrector.slurm NP1free_R1_001.fastq.gz NP1free_R2_001.fastq.gz
```

## Step 3: Discard read pairs with one damaged read
# download python script here: https://github.com/harvardinformatics/TranscriptomeAssemblyTools

#!/bin/bash
#SBATCH -J filter_reads_$1           # Job name
#SBATCH -n 1                         # Use 1 core for the job
#SBATCH -t 06:00:00                  # Runtime in HH:MM:SS
#SBATCH -p shared            # Partition to submit to
#SBATCH --mem=20000                   # Memory per node in MB
#SBATCH -o rmunfix_%A.o              # File to which STDOUT will be written
#SBATCH -e rmunfix_%A.e              # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=pvaelli@fas.harvard.edu # Email to which notifications will be sent

module purge
module load python/2.7.8-fasrc01

#path to directory
cwd=$(pwd)

python $cwd/FilterUncorrectablePEfastq.py -1 $1 -2 $2 -s $3


# Submit the job with following line:
sbatch filter_reads.sh Body_R1_001.cor.fq.gz Body_R2_001.cor.fq.gz Body_summary

# Note: output files will have prefix unfixrm_ and will be unzipped


## Step 4: Trim seqs with Trimgalore 

# note: this script uses an older version of Trimgalore, which is currently in my folder in the shared bellono_lab directory
# trimgalore must be in same directory with files. This job ran for ~15 hrs

# submission script:

#!/bin/bash 
#SBATCH -n 1 #Number of cores
#SBATCH -t 24:00:00 #Runtime in minutes
#SBATCH -p shared #Partition to submit to
#SBATCH -e trimgalore_%A.e
#SBATCH -o trimgalore_%A.o
#SBATCH --mem=3000 #Memory per node in MB
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=<pvaelli@fas.harvard.edu>  # Email to which notifications will be sent

echo -n "Starting job on "
date

module load cutadapt/1.8.1-fasrc01
cwd=$(pwd)

for file in *1_001.cor.fq; do
  out=${file%1_001.cor.fq}
  
$cwd/trim_galore --paired --phred33 --length 36 -q 0 --stringency 1 -e 0.1 ./${out}1_001.cor.fq ./${out}2_001.cor.fq
done


## Step 5: Map reads to SILVA rRNA database to remove contamination
# this method could be used to map reads to any black list (parasites, virus, fungi, etc)

# Make a bowtie2 index from SILVA rRNA LSU and SSU parc files. 
# concatenate files:
cat file1 file2

# convert U to T (RNA to DNA seq):
sed '/^[^>]/s/U/T/g' file.fasta >newfile.fasta

# Submit script to run bowtie2 and make index

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12 #Number of cores
#SBATCH -t 12:00:00  #Runtime in minutes
#SBATCH -p shared  #Partition to submit to
#SBATCH --mem=48000  #Memory per node in MB
#SBATCH -e bowtie2_%A.e
#SBATCH -o bowtie2_%A.o
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=<pvaelli@fas.harvard.edu>  # Email to which notifications will be sent


module purge
module load bowtie2/2.3.2-fasrc02

bowtie2-build --large-index --threads 12 -f SILVA_rRNA_clean.fasta SILVA_rRNA


# After the index is finished, include the path to the directory with all index files:

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12 #Number of cores
#SBATCH -t 12:00:00  #Runtime in minutes
#SBATCH -p shared  #Partition to submit to
#SBATCH --mem=48000  #Memory per node in MB
#SBATCH -e silvabt2_%A.e
#SBATCH -o silvabt2_%A.o
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=<pvaelli@fas.harvard.edu>  # Email to which notifications will be sent

# $1 = full path to your silva database; do not include the fasta suffix (fa,fasta,etc.)
# $2 = R1 fastq file
# $3 = R2 fastq file
# $4 = sample_id (no spaces)

singularity exec --cleanenv /n/singularity_images/informatics/trinityrnaseq/trinityrnaseq.v2.11.0.simg bowtie2 --quiet --very-sensitive-local --phred33  -x $1 -1 $2 -2 $3 --threads 12 --met-file ${4}_bowtie2_metrics.txt --al-conc-gz blacklist_paired_aligned_${4}.fq.gz --un-conc-gz blacklist_paired_unaligned_${4}.fq.gz  --al-gz blacklist_unpaired_aligned_${4}.fq.gz --un-gz blacklist_unpaired_unaligned_${4}.fq.gz


sbatch rRNA_cleanup.sh /n/holyscratch01/bellono_lab/Users/pvaelli/plastid_cells/rRNA_database/SILVA_rRNA P1poly_R1_trimmed.fq.gz P1poly_R2_trimmed.fq.gz P1_poly
sbatch rRNA_filtering.slurm /n/holyscratch01/bellono_lab/Users/pvaelli/penicillus/rRNA_database/SILVA_rRNA Algae_1.cor.trimmed.fq.gz Algae_2.cor.trimmed.fq.gz Algae_filtered
sbatch rRNA_filter.sh /n/holyscratch01/bellono_lab/Users/pvaelli/Eclarki/rRNA_database/SILVA_rRNA unfixrm_Parapodia_R1_001.cor_val_1.fq unfixrm_Parapodia_R2_001.cor_val_2.fq Parapodia_filtered

sbatch rRNA_cleanup.sh /n/holyscratch01/bellono_lab/Users/pvaelli/rRNA_database/SILVA_rRNA

# read pairs for which neither read mapped to the rRNA database,i.e. "paired_unaligned" reads, specified after the --un-conc-gz flag.


## Step 6: Quality check and remove overrepresented seqs

# run FASTQC:

# submission script:

#!/bin/bash
#SBATCH -p shared       # Partition to submit to
#SBATCH -n 1                   # Number of cores
#SBATCH -t 1-8:00               # Runtime in days-hours:minutes
#SBATCH --mem 2000              # Memory in MB
#SBATCH -J FastQC               # job name
#SBATCH -o FastQC.%A.out        # File to which standard out will be written
#SBATCH -e FastQC.%A.err        # File to which standard err will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=<pvaelli@fas.harvard.edu>  # Email to which notifications will be sent

source new-modules.sh
module purge
module load fastqc/0.11.5-fasrc01

for file in *.fq.gz; do
        fastqc $file
done

# Unzip each new directory and rename/copy the txt file up to the fastq files
# Then run this python script which reads these txt files and removes overrep kmers from the fastq files

#!/bin/bash
#SBATCH -J remove_overrep_$1                # Job name
#SBATCH -n 1                         # Use 1 core for the job
#SBATCH -t 04:00:00                  # Runtime in HH:MM:SS
#SBATCH -p shared            # Partition to submit to
#SBATCH --mem=5000                   # Memory per node in MB
#SBATCH -o remove_overrep_%A.o              # File to which STDOUT will be written
#SBATCH -e remove_overrep_%A.e              # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=pvaelli@fas.harvard.edu # Email to which notifications will be sent

module purge
module load python/2.7.8-fasrc01

python /n/holyscratch01/bellono_lab/Users/pvaelli/plastid_cells/filtered_reads/RemoveFastqcOverrepSequenceReads.py -1 $1 -2 $2 -fql $3 -fqr $4

# sbatch overrep_seq.sh NP1poly_filtered_R1.fq.gz NP1poly_filtered_R2.fq.gz NP1poly_filtered_R1_fastqc.txt NP1poly_filtered_R2_fastqc.txt



## Step 7: Run Trinity


# Step 7A Harvard FAS method:

sbatch trinity.sh --left NP1free_R1_001_val_1.fq.gz,NP1poly_R1_001_val_1.fq.gz,P1free_R1_001_val_1.fq.gz,P1poly_R1_001_val_1.fq.gz --right NP1free_R2_001_val_2.fq.gz,NP1poly_R2_001_val_2.fq.gz,P1free_R2_001_val_2.fq.gz,P1poly_R2_001_val_2.fq.gz
sbatch trinity.sh --left blacklist_paired_unaligned_Algae_filtered.fq.1.gz --right blacklist_paired_unaligned_Algae_filtered.fq.2.gz
sbatch trinity.sh --left Body_final_R1.fq.gz,Parapodia_final_R1.fq.gz --right Body_final_R2.fq.gz,Parapodia_final_R2.fq.gz


#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
#SBATCH --partition=shared
#SBATCH -e trinity_%A.err            # File to which STDERR will be written
#SBATCH -o trinity_%A.out           # File to which STDOUT will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=pvaelli@fas.harvard.edu # Email to send notifications
set -o nounset -o errexit -o xtrace

########################################
# parameters
########################################
readonly SINGULARITY_IMAGE=/n/singularity_images/informatics/trinityrnaseq/trinityrnaseq.v2.11.0.simg
readonly TRINITY_OUT_DIR=trinity_out_dir
# To see all options:
#    singularity exec --cleanenv ${SINGULARITY_IMAGE}  Trinity --show_full_usage_info
readonly TRINITY_OPTIONS="--output ${TRINITY_OUT_DIR} --max_memory $((8*$(ulimit -m)/(1024**2)/10))G --CPU ${SLURM_CPUS_ON_NODE} --no_normalize_reads --seqType fq $@"

########################################
# ... don't modify below here ...
if [ ! -s "${TRINITY_OUT_DIR}/read_partitions.img" ]
then
  mkdir -p "${TRINITY_OUT_DIR}"
  readonly tmpdir=$(mktemp -d)
  mkdir -m 777 -p ${tmpdir}/upper ${tmpdir}/work
  truncate -s 2T "${TRINITY_OUT_DIR}/read_partitions.img"
  singularity exec --cleanenv ${SINGULARITY_IMAGE} mkfs.ext3 -d "${tmpdir}" "${TRINITY_OUT_DIR}/read_partitions.img"
  singularity exec --cleanenv --overlay ${TRINITY_OUT_DIR}/read_partitions.img ${SINGULARITY_IMAGE} mkdir /read_partitions
  ln -sf /read_partitions ${TRINITY_OUT_DIR}/read_partitions
  rm -rf "${tmpdir}"
fi

# if on a bigmem node, stop after inchworm
case ${SLURM_JOB_PARTITION} in
  bigmem) no_run_chrysalis='--no_run_chrysalis' ;;
       *) no_run_chrysalis='' ;;
esac

srun -n 1 env time -v singularity exec \
                           --cleanenv \
                           --no-home \
                           --overlay ${TRINITY_OUT_DIR}/read_partitions.img \
                           "${SINGULARITY_IMAGE}" \
  Trinity ${TRINITY_OPTIONS} ${no_run_chrysalis}
  
  

 
# Step 7B: Bigmem method:

# Subsample reads to reduce redundancy and computational burden
# My four samples are between 70-100 million reads each, aiming for 100 million reads total, so subsample to 25 million each
# Note: this command will produce uncompressed files, even though extension is still fastq.gz
# You need to remove .gz from each file name and re-compress using gzip command or Trinity will not work

# submission script:

#!/bin/bash 
#SBATCH -n 1 #Number of cores
#SBATCH -t 24:00:00 #Runtime in minutes
#SBATCH -p shared #Partition to submit to
#SBATCH -e downsample_%A.err
#SBATCH -o downsample_%A.out
#SBATCH --mem=10000 #Memory per node in MB
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=<pvaelli@fas.harvard.edu>  # Email to which notifications will be sent

module purge
module load seqtk/1.2-fasrc01

for file in *.fq.gz; do
	  seqtk sample -s100 $file 25000000 > sub/$file
done

# Check that it worked. Count the number of reads in original vs. subsampled files

zipped file: echo $(zcat yourfile.fq.gz|wc -l)/4|bc
unzipped: echo $(cat yourfile.fq|wc -l)/4|bc

# example
echo $(zcat P1poly_R1_001_val_1.fq.gz|wc -l)/4|bc


# Bigmem sub script method:

#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p bigmem                   # may consider running on a bigmem node for large dataset
#SBATCH -e trinity_%A.err            # File to which STDERR will be written
#SBATCH -o trinity_%A.out           # File to which STDOUT will be written
#SBATCH -J trinity_%A               # Job name
#SBATCH --mem=225000                # Memory requested
#SBATCH --time=10-00:00:00              # Runtime in D-HH:MM:SS
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=pvaelli@fas.harvard.edu # Email to send notifications to

module purge
module load trinityrnaseq/2.4.0-fasrc02
# $1 = comma-separated list of R1 files
# $2 = comma-separated list of R2 files
# $3 = name of output directory Trinity will create to store results. This must include Trinity in the name, otherwise the job will terminate

Trinity --seqType fq --max_memory 200G --min_kmer_cov 1 --CPU 24 --left blacklist_paired_unaligned_Algae_filtered.fq.1.gz --right blacklist_paired_unaligned_Algae_filtered.fq.2.gz --output /n/holyscratch01/bellono_lab/Users/pvaelli/penicillus/trinity_out



# Step 8: Transdecoder ORF finder

#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 12000
#SBATCH -p shared
#SBATCH -J transdecoder               # job name
#SBATCH -o transdecoder.out
#SBATCH -e transdecoder.err
#SBATCH -t 0-24:00               # Runtime in days-hours:minutes
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=<pvaelli@fas.harvard.edu>  # Email to which notifications will be send to
 
source new-modules.sh
module load TransDecoder/5.3.0-fasrc01
 
TRANS_DATA=$(pwd)

TransDecoder.LongOrfs -t $TRANS_DATA/*Trinity.fasta
TransDecoder.Predict -t $TRANS_DATA/*Trinity.fasta --single_best_only

# Step 9: Swissprot annotation using Diamond

# Need to first make a diamond index with FASTA reference database:

#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p shared
#SBATCH -e diamond_%A.e
#SBATCH -o diamond_%A.o
#SBATCH -J diamond
#SBATCH --mem=12000
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pvaelli@fas.harvard.edu


module purge
module load gcc/9.3.0-fasrc01
module load diamond/2.0.4-fasrc01

diamond makedb --in swissprot.fasta -d swissprot_diamond

# Then use database for Diamond:

#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p shared
#SBATCH -e diamond_%A.e
#SBATCH -o diamond_%A.o
#SBATCH -J diamond
#SBATCH --mem=12000
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pvaelli@fas.harvard.edu


module purge
module load gcc/9.3.0-fasrc01
module load diamond/2.0.4-fasrc01

DATA=*.Trinity.fasta.transdecoder.pep
DATABASE=/n/holyscratch01/bellono_lab/Users/pvaelli/swissprot_diamond
OUT=$(echo ${DATA} | sed 's/\..*//')

# set -k 1 flag to retrieve one target seq per query; otherwise default=25

diamond blastp -q $DATA -d $DATABASE -o $OUT.diamond.tsv -k 1 -f 6 qseqid stitle pident evalue length mismatch gapopen qstart qend sstart send bitscore --ultra-sensitive



# Step 10: Gene quantification using Kallisto read mapping

# first make kallisto index with cds file from transdecoder:

#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p shared
#SBATCH -e kallisto_index_%A.e
#SBATCH -o kallisto_index_%A.o
#SBATCH -J kallisto_index
#SBATCH --mem=12000
#SBATCH --time=12:00:00

module purge
module load kallisto/0.45.1-fasrc01

# make kallisto index
# replace index with an index name
# replace transcriptome.fa with Trinity assembly
# kallisto index -i index transcriptome.fa

DATA=*.Trinity.fasta.transdecoder.cds
OUT=$(echo ${DATA} | sed 's/\..*//')


kallisto index -i $OUT.kallisto $DATA


# Then submit each pair of reads:

#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p shared
#SBATCH -e kallisto%A.e
#SBATCH -o kallisto%A.o
#SBATCH -J kallisto
#SBATCH --mem=12000
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pvaelli@fas.harvard.edu


module purge
module load kallisto/0.45.1-fasrc01

# replace index with appropriate index name from above
# replace mapping/1 with an output folder (diff one per thing to map)
# replace the trimmed/fastq.gz files with paired R1 and R2 of input files (trimmed_
# kallisto quant -i index -o mapping/1 trimmed/1_1.R1.fastq.gz trimmed/1_2.R2.fastq.gz

DATABASE=*.kallisto

for file in *1.fq.gz; do
  out=${file%*1.fq.gz}
  
NAME=$(echo ${out} | sed 's/\_.*//')
kallisto quant -i $DATABASE -o $NAME ${out}1.fq.gz ${out}2.fq.gz

done



