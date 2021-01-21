# SKSV
ultrafast structural variation detection from circular consensus sequencing reads

## Getting Start
    git clone --recursive https://github.com/ydLiu-HIT/SKSV.git
    cd SKSV/skeleton
    make   ## built deBGA for RdBG-index
    cd ..
    make   ## built deSALT for alignment
    
    ./deSALT index ref.fa index_route
    ./deSALT aln index_route read.fq
    
    or 
    run deSALT directly in the same folder (Executable programs have been built in advance.)
    
    ## install by conda
    conda install -c bioconda desalt
