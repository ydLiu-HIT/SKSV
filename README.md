# SKSV
ultrafast structural variation detection from circular consensus sequencing reads

## Getting Start
    git clone https://github.com/ydLiu-HIT/SKSV.git
    cd SKSV/skeleton/
    make   ## built ./SKSV-skeleton
    cd ..
    
    ./SKSV index ref.fa index_route   # build index 
    ./SKSV aln index_route read.fq    # skeleton-alignment
    ./SKSV call in.svseg ref.fa out.vcf work_dir # call variants using in.svseg file
        
    ## add to enviroment path
    vim ~/.bashrc
    export PATH="path:$PATH"  # path is the folder where you download SKSV
    
       
## Introduction

Circular consensus sequencing (CCS) reads are promising for the comprehensive detection of structural variants (SVs). However, alignment-based SV calling pipelines are computationally intensive due to the generation of complete read-alignments and its post-processing. Herein, we propose a **SK**eleton-based analysis toolkit for **S**tructural **V**ariation detection (SKSV).  Benchmarks on real and simulated datasets demonstrate that SKSV has an order of magnitude of faster speed than state-of-the-art SV calling ap-proaches, moreover, it enables to well-handle various types of SVs with higher F1 scores.



## Synopsis
Reference genome indexing
```
SKSV index ref.fa <index_route>
```
	
Skeleton alignment
```
SKSV aln <index_route> read.fa/fq -o skeleton.svseg
```

Variant calling
```
SKSV call skeleton.svseg ref.fa output.vcf workdir
```

## Commands and options
**Genome index**
```
Usage:     SKSV-skeleton index [Options] <Reference> <HashIndexDir>
<Reference>            Sequence of reference genome, in FASTA format
<HashIndexDir>         The directory to store deBGA index
``` 

**Skeleton alignment**
```
Usage:     SKSV-skeleton aln [options] <index_route> <read.fa/fq>

    <index_route>                 The path of RdBG index.
    <read.fq/fa>                  The input reads in fasta or fastq format.

Algorithm options:

    -t --thread           [INT]    Number of threads. [4]
    -k --index_kmer       [INT]    RdBG-index kmer length. [22]
    -s --seeding_kmer     [INT]    K-mer length of seeding process (no long than RdBG-index kmer). [22]
    -d --hash_kmer        [INT]    K-mer length of local hash process. [7]
    -B --batch_size       [INT]    The number of reads to be processed in one loop. [100000]
    -l --sv_lim           [INT]    The minimum length to consider a structure variant. [30]
    -m --error_model      [INT]    Data type for different error model. [1]
    -w --block            [INT]    The window length when extend node in the skeleton. [30]
    -x --edit_dis         [INT]    The maximum edit distance for Landua-Vishkin. [0]
    -e --seed_step        [INT]    The interval of seeding. [20]
    -M --extend_mx        [INT]    The maximum length for an anchor. [5000]
    -N --top_N            [INT]    Max allowed skeleton paths. [4]
    -g --rd_mx_gap        [INT]    Maximum allowed distance in read that two anchor can be connected. [1000]
    -G --rf_mx_gap        [INT]    Maximum allowed distance in reference that two anchor can be connected. [1000]
    -Y --min_chain_score  [INT]    The minimum skeleton chaining score to be processed. [100]
    -p --secondary_ratio  [FLAOT]  Minimum secondary-to-primary alignment score ratio. [0.8]
    -P --overlap_ratio    [FLAOT]  Minimum overlap ratio to consider tow skeletons one as primary one as secondary. [0.4]

Output options:

    -o --output    	  [STR]    Output file for SV signatures in svseg format. [sk.svseg]
    -h --help                      Show detailed usage.
```

**Variants calling**
```
positional arguments:
  input                 Skeletons from SKSV-skeleton in svseg fromat.
  ref                   The reference genome in fasta format.
  output                Output VCF format file.
  work_dir              Work-directory for distributed jobs

optional arguments:
  -h, --help            show this help message and exit
  --version, -v         show program's version number and exit
  -t THREADS, --threads THREADS
                        Number of threads to use.[16]
  -b BATCHES, --batches BATCHES
                        Batch of genome segmentation interval.[10000000]
  -S SAMPLE, --sample SAMPLE
                        Sample name/id
  --retain_work_dir     Enable to retain temporary folder and files.

Collection of SV signatures:
  -p MAX_SPLIT_PARTS, --max_split_parts MAX_SPLIT_PARTS
                        Maximum number of split segments a read may be aligned before it is ignored.[7]
  --merge_del_threshold MERGE_DEL_THRESHOLD
                        Maximum distance of deletion signals to be merged.[500]
  --merge_ins_threshold MERGE_INS_THRESHOLD
                        Maximum distance of insertion signals to be merged.[500]

Generation of SV clusters:
  -s MIN_SUPPORT, --min_support MIN_SUPPORT
                        Minimum number of reads that support a SV to be reported.[10]
  -l MIN_SIZE, --min_size MIN_SIZE
                        Minimum size of SV to be reported.[30]
  -L MAX_SIZE, --max_size MAX_SIZE
                        Maximum size of SV to be reported.[100000]
  -sl MIN_SIGLENGTH, --min_siglength MIN_SIGLENGTH
                        Minimum length of SV signal to be extracted.[10]

Computing genotypes:
  --genotype            Enable to generate genotypes.
  --gt_round GT_ROUND   Maximum round of iteration for alignments searching if perform genotyping.[500]
  
Generate allele sequence:
  --read READ           Raw reads in bgzip fasta/fastq format with index by samtools for extract inserted sequence for insertions.
  --print_allele_seq    Enable to print inserted and deleted sequence for an insertion and deletion variants. If --print_allele_seq was set, --read is needed.

Advanced:
  --max_cluster_bias_INS MAX_CLUSTER_BIAS_INS
                        Maximum distance to cluster read together for insertion.[800]
  --diff_ratio_merging_INS DIFF_RATIO_MERGING_INS
                        Do not merge breakpoints with basepair identity more than [0.9] for insertion.
  --max_cluster_bias_DEL MAX_CLUSTER_BIAS_DEL
                        Maximum distance to cluster read together for deletion.[1000]
  --diff_ratio_merging_DEL DIFF_RATIO_MERGING_DEL
                        Do not merge breakpoints with basepair identity more than [0.5] for deletion.
  --max_cluster_bias_INV MAX_CLUSTER_BIAS_INV
                        Maximum distance to cluster read together for inversion.[500]
  --max_cluster_bias_DUP MAX_CLUSTER_BIAS_DUP
                        Maximum distance to cluster read together for duplication.[500]
  --max_cluster_bias_TRA MAX_CLUSTER_BIAS_TRA
                        Maximum distance to cluster read together for translocation.[50]
  --diff_ratio_filtering_TRA DIFF_RATIO_FILTERING_TRA
                        Filter breakpoints with basepair identity less than [0.6] for translocation.
```

## Datasets 
We implemented the benchmarks on a server with AMD Ryzen 3950X CPU and 120GB RAM, running Linux Ubuntu 16.04.

The real HG002 PacBio CCS reads were downloaded from: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_SequelII_CCS_11kb/reads/. 

The GIAB ground truth set and corresponding high confidence region set were downloaded from: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/.

The new draft difficult medically relevant gene Structural Variant benchmark from GIAB (Zook, et al., 2016) that includes some SVs in segmental duplications and corresponding benchmark bed files were download from https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/

Truvari (v2.0.0, https://github.com/spiralgenetics/truvari) was used to assess the precision, recall, and F1 score (Zook, et al., 2020)..

## Reference

Zook, J.M., et al. Extensive sequencing of seven human genomes to characterize benchmark reference materials. Sci Data 2016;3.

Zook, J.M., et al. A robust benchmark for detection of germline large deletions and insertions. Nat Biotechnol 2020.


## Contact
For advising, bug reporting and requiring help, please post on **[Github Issues](https://github.com/ydLiu-HIT/SKSV/issues)** or contact ydliu@hit.edu.cn.

