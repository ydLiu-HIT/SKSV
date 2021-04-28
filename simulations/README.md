Due to the lack of ground truth for inversions, duplications and translocations. We implemented a simulation benchmark similar to that of cuteSV study. 44 inversion were extracted from CHM1 samples callsets (Huddleston, et al., 2017) (nstd137 in dbVAR database) and 3712 and 380 non-overlapping duplications and translocations were extracted from KWS1 sample callsets (Alsmadi, et al., 2014) (nstd137 in dbVAR database). Then the three types of SVs were respectively integrated into humna reference genome (hs37d5) to build three in silico donor genome to generate simulated datasets using VISOR simulator.

The extracted SV events are transformed into BED format according to the requirement of VISOR HACk, i.e. sim_inv.bed, sim_tra.bed, sim_dup.bed. LASeR.bed describes the genomic region and depth to simulate.

Moreover, another 518 no-overlapping deletions in chromosome 2 (sim_del_chr2.bed) was extracted from CHM1 samples callsets to generate a number of simulated datasets in 1%-15% error rates, to assess the ability of SKSV on various sequencing error rates.

### Generate a simulated reference with haplotype-specific structural variants

```
for type in {inv,dup,tra,}; do VISOR HACk -g hs37d5.fasta -bed sim_${type}.bed -o donor_genome_${type}; done
```

### Generate simulated reads using the simulated reference genome
```
for type in {inv,dup,tra,}; do VISOR LASeR -g hs37d5.fa -s donor_genome_${type}/ -bed LASeR.bed -o data_${type}_30x -c 30 --threads 16 --noaddtag -a 0.99 -l 12000 -r 10:60:30 --readstype PB --ccs; done

```
