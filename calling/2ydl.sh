BINDIR=$1

INSIGS=$2
FOLDER=$3
VCFNAME=$4
THREAD=$5
GENOTYPE=$6

vcf=${FOLDER}"/"${VCFNAME}".vcf"
vcfgz=${vcf}".gz"

pypy3 ${BINDIR}/load_rst-v1.6.py ${INSIGS} ${FOLDER}/signals.txt
cat ${FOLDER}/signals.txt | grep DEL | sort -u | sort -k 2,2 -k 3,3n > ${FOLDER}/DEL.sigs
cat ${FOLDER}/signals.txt | grep INS | sort -u | sort -k 2,2 -k 3,3n > ${FOLDER}/INS.sigs
#python3 /home/ydliu/deVar/sv_bench/cuteSV-devar-v1.6/cuteSV-deVar.py /home/ydliu/MLSV/data/HG002.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam  ${vcf} ./ -s 3 -l 30 --max_cluster_bias_INS 500 --diff_ratio_merging_INS 0.65 --diff_ratio_filtering_INS 0.65 --diff_ratio_filtering_DEL 0.35 -mi 500 -md 500

##for ccs
#python ${BINDIR}/cuteSV-deVar.py /data/0/ydliu/HG002/HG002.SequelII.11kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam  ${vcf} ${FOLDER} -s 3 -l 30 --max_cluster_bias_INS 800 --diff_ratio_merging_INS 0.65 --diff_ratio_filtering_INS 0.65 --diff_ratio_filtering_DEL 0.35 -t ${THREAD}

python ${BINDIR}/cuteSV-deVar.py ${INSIGS} /data/0/ydliu/HG002/HG002.SequelII.11kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam  /data/0/ydliu/MLSV/REF/hs37d5.fa ${vcf} ${FOLDER} -s 3 -l 30 --max_cluster_bias_INS 800 --diff_ratio_merging_INS 0.5 --diff_ratio_filtering_INS 0.5 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.4 --diff_ratio_filtering_DEL 0.5 -t ${THREAD} ${GENOTYPE}

#for clr
#python ${BINDIR}/cuteSV-deVar.py /data/0/ydliu/HG002/HG002.SequelII.11kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam  ${vcf} ${FOLDER} -s 8 -l 30 --max_cluster_bias_INS 800 --max_cluster_bias_DEL 300 --diff_ratio_merging_INS 0.2 --diff_ratio_filtering_INS 0.6 --diff_ratio_filtering_DEL 0.7 -t ${THREAD}
#python ${BINDIR}/cuteSV-deVar.py ${INSIGS} /data/0/ydliu/HG002/HG002.SequelII.11kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam /data/0/ydliu/MLSV/REF/hs37d5.fa ${vcf} ${FOLDER} -s 8 -l 30 --max_cluster_bias_INS 800 --diff_ratio_merging_INS 0.5 --diff_ratio_filtering_INS 0.8 --max_cluster_bias_DEL 300 --diff_ratio_merging_DEL 0.2 --diff_ratio_filtering_DEL 0.6 -t ${THREAD}

echo ${vcf}
echo ${vcfgz}
bgzip -c ${vcf} > ${vcfgz}
tabix ${vcfgz}

