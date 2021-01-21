# SKSV
ultrafast structural variation detection from circular consensus sequencing reads

## Getting Start
    git clone --recursive https://github.com/ydLiu-HIT/SKSV.git
    cd SKSV/skeleton/
    make   ## built ./SKSV-skeleton
    cd ..
    
    ./SKSV index ref.fa index_route   # build index 
    ./SKSV aln index_route read.fq    # skeleton-alignment
    ./SKSV call in.svseg ref.fa out.vcf word_dir # call variants using in.svseg file
        
    ## add to enviroment path
    vim ~/.bashrc
    export PATH="path:$PATH"  # path is the folder where you download SKSV
    
       
## Introduction
The rising circular consensus sequencing (CCS) long-reads promise unprecedentedly comprehensive structural variants detection. However, the state-of-the-art alignment-based SV calling pipelines are computationally intensive due to the generation of complete read-alignments and its post-processing. Herein, we propose a **Sk**eleton-based analysis toolkit for **S**tructural **V**ariation detection (SKSV). Compared with the alignment-based SV calling pipeline on a HG002 real CCS dataset, SKSV achieves better F1 scores for discovering insertions and deletions with an order of magnitude of acceleration.



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
**index**
Usage:     SKSV-skeleton index [Options] <Reference> <HashIndexDir>
<Reference>            Sequence of reference genome, in FASTA format
<HashIndexDir>         The directory to store deBGA index
    

**aln**
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
    -x --edit_dis         [INT]    The maximum edit distance for Landua-Vishkin. [3]
    -e --seed_step        [INT]    The interval of seeding. [20]
    -M --extend_mx        [INT]    The maximum length for an anchor. [200]
    -N --top_N            [INT]    Max allowed skeleton paths. [3]
    -g --rd_mx_gap        [INT]    Maximum allowed distance in read that two anchor can be connected. [1500]
    -G --rf_mx_gap        [INT]    Maximum allowed distance in reference that two anchor can be connected. [2000]
    -Y --min_chain_score  [INT]    The minimum skeleton chaining score to be processed. [100]
    -p --secondary_ratio  [FLAOT]  Minimum secondary-to-primary alignment score ratio. [0.9]
    -P --overlap_ratio    [FLAOT]  Minimum overlap ratio to consider tow skeletons one as primary one as secondary. [0.4]

Output options:

    -o --output    [STR]    Output file for SV signatures in svseg format. [sk.svseg]
    -h --help                      Show detailed usage.



