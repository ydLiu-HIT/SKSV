# SKSV
ultrafast structural variation detection from circular consensus sequencing reads

## Getting Start
    git clone https://github.com/ydLiu-HIT/SKSV.git
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
  -q MIN_MAPQ, --min_mapq MIN_MAPQ
                        Minimum mapping quality value of alignment to be taken into account.[20]
  -r MIN_READ_LEN, --min_read_len MIN_READ_LEN
                        Ignores reads that only report alignments with not longer than bp.[500]
  -md MERGE_DEL_THRESHOLD, --merge_del_threshold MERGE_DEL_THRESHOLD
                        Maximum distance of deletion signals to be merged.[0]
  -mi MERGE_INS_THRESHOLD, --merge_ins_threshold MERGE_INS_THRESHOLD
                        Maximum distance of insertion signals to be merged.[100]

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

Advanced:
  --max_cluster_bias_INS MAX_CLUSTER_BIAS_INS
                        Maximum distance to cluster read together for insertion.[100]
  --diff_ratio_merging_INS DIFF_RATIO_MERGING_INS
                        Do not merge breakpoints with basepair identity more than [0.2] for insertion.
  --diff_ratio_filtering_INS DIFF_RATIO_FILTERING_INS
                        Filter breakpoints with basepair identity less than [0.6] for insertion.
  --max_cluster_bias_DEL MAX_CLUSTER_BIAS_DEL
                        Maximum distance to cluster read together for deletion.[200]
  --diff_ratio_merging_DEL DIFF_RATIO_MERGING_DEL
                        Do not merge breakpoints with basepair identity more than [0.3] for deletion.
  --diff_ratio_filtering_DEL DIFF_RATIO_FILTERING_DEL
                        Filter breakpoints with basepair identity less than [0.7] for deletion.
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
We implemented the benchmarks on a server with Intel Xeon E4280 2.0GHZ CPU and 1 Terabytes RAM, running Linux Ubuntu 16.04. 

The real HG002 PacBio CCS reads were downloaded from: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_SequelII_CCS_11kb/reads/. 

The GIAB ground truth set and corresponding high confidence region set were downloaded from: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/.

Truvari (v2.0.0, https://github.com/spiralgenetics/truvari) was used to assess the precision, recall, and F1 score.

## Contact
For advising, bug reporting and requiring help, please post on **[Github Issues](https://github.com/ydLiu-HIT/SKSV/issues)** or contact ydliu@hit.edu.cn.

