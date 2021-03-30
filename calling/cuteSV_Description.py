''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  cuteSV_Description.py
 * @author: tjiang
 * @date: Apr 26th 2020
 * @version V1.0.10   
'''
import argparse

VERSION_1 = '1.0.2'
VERSION = '1.0.10'


class cuteSVdp(object):
    '''
    Detailed descriptions of cuteSV version and its parameters.

    '''

    USAGE = """\
	Long read based fast and accurate SV detection with cuteSV.
	
	Current version: v%s
	Author: Tao Jiang, Yadong Liu
	Contact: tjiang@hit.edu.cn


	Suggestions:

	For PacBio CCS(HIFI) data:
		--max_cluster_bias_INS	        800	
		--diff_ratio_merging_INS	0.9
                --max_cluster_bias_DEL          1000
                --diff_ratio_merging_DEL        0.5
                --merge_del_threshold           500
                --merge_ins_threshold           500


	""" % (VERSION)

    # MinSizeDel = 'For current version of cuteSV, it can detect deletions larger than this size.'


def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="cuteSV",
                                     description=cuteSVdp.USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--version', '-v',
                        action='version',
                        version='%(prog)s {version}'.format(version=VERSION))

    # **************Parameters of input******************
    parser.add_argument("input",
                        type=str,
                        help="Skeletons from SKSV-skeleton in svseg fromat.") 
    parser.add_argument("ref",
                        type=str,
                        help="The reference genome in fasta format.")
    parser.add_argument('output',
                        type=str,
                        help="Output VCF format file.")
    parser.add_argument('work_dir',
                        type=str,
                        help="Work-directory for distributed jobs")

    # ************** Other Parameters****************** 
    parser.add_argument('-t', '--threads',
                        help="Number of threads to use.[%(default)s]",
                        default=16,
                        type=int)
    parser.add_argument('-b', '--batches',
                        help="Batch of genome segmentation interval.[%(default)s]",
                        default=10000000,
                        type=int)
    # The description of batches needs to improve.
    parser.add_argument('-S', '--sample',
                        help="Sample name/id",
                        default="NULL",
                        type=str)

    parser.add_argument('--retain_work_dir',
                        help="Enable to retain temporary folder and files.",
                        action="store_true")

    parser.add_argument('--report_readid',
                        help="Enable to report supporting read ids for each SV.",
                        action="store_true")

    # **************Parameters in signatures collection******************
    GroupSignaturesCollect = parser.add_argument_group(
        'Collection of SV signatures')
    GroupSignaturesCollect.add_argument('-r', '--min_read_len',
                                        help="Ignores reads that only report alignments with not longer than bp.[%(default)s]",
                                        default=500,
                                        type=int)
    GroupSignaturesCollect.add_argument('--merge_del_threshold',
                                        help="Maximum distance of deletion signals to be merged. In our paper, I used -md 500 to process HG002 real human sample data.[%(default)s]",
                                        default=0,
                                        type=int)
    GroupSignaturesCollect.add_argument('--merge_ins_threshold',
                                        help="Maximum distance of insertion signals to be merged. In our paper, I used -mi 500 to process HG002 real human sample data.[%(default)s]",
                                        default=100,
                                        type=int)
    # The min_read_len in last version is 2000.
    # signatures with overlap need to be filtered

    # **************Parameters in clustering******************
    GroupSVCluster = parser.add_argument_group('Generation of SV clusters')
    GroupSVCluster.add_argument('-s', '--min_support',
                                help="Minimum number of reads that support a SV to be reported.[%(default)s]",
                                default=10,
                                type=int)
    GroupSVCluster.add_argument('-l', '--min_size',
                                help="Minimum size of SV to be reported.[%(default)s]",
                                default=30,
                                type=int)
    GroupSVCluster.add_argument('-L', '--max_size',
                                help="Maximum size of SV to be reported.[%(default)s]",
                                default=100000,
                                type=int)
    GroupSVCluster.add_argument('-sl', '--min_siglength',
                                help="Minimum length of SV signal to be extracted.[%(default)s]",
                                default=10,
                                type=int)

    # **************Parameters in genotyping******************
    GroupGenotype = parser.add_argument_group('Computing genotypes')
    GroupGenotype.add_argument('--genotype',
                               help="Enable to generate genotypes.",
                               action="store_true")
    GroupGenotype.add_argument('--gt_round',
                               help="Maximum round of iteration for alignments searching if perform genotyping.[%(default)s]",
                               default=500,
                               type=int)

    # **************Parameters in generating allele sequence******************
    GroupSequence = parser.add_argument_group('Generate allele sequence')
    GroupSequence.add_argument("--read",
                        type=str,
                        default="None",
                        help="Raw reads in bgzip fasta/fastq format with index by samtools for extract inserted sequence for insertions.")
    GroupSequence.add_argument('--print_allele_seq',
                               help="Enable to print inserted/deleted sequence for an insertion/deletion variant. If --print_allele_seq was set, --read is needed.",
                               action="store_true")

    # GroupGenotype.add_argument('--hom',
    # 	help = "Threshold on allele frequency for homozygous.[%(default)s]",
    # 	default = 0.8,
    # 	type = float)
    # GroupGenotype.add_argument('--het',
    # 	help = "Threshold on allele frequency for heterozygous.[%(default)s].",
    # 	default = 0.2,
    # 	type = float)

    # Just a parameter for debug.
    # Will be removed in future.
    # GroupSVCluster.add_argument('--preset',
    # 	help = "Parameter presets for different sequencing technologies (pbclr/pbccs/ont).[%(default)s]",
    # 	default = "pbccs",
    # 	type = str)

    # **************Advanced Parameters******************
    GroupAdvanced = parser.add_argument_group('Advanced')

    # ++++++INS++++++
    GroupAdvanced.add_argument('--max_cluster_bias_INS',
                               help="Maximum distance to cluster read together for insertion.[%(default)s]",
                               default=800,
                               type=int)
    GroupAdvanced.add_argument('--diff_ratio_merging_INS',
                               help="Do not merge breakpoints with basepair identity more than [%(default)s] for insertion.",
                               default=0.9,
                               type=float)
    # ++++++DEL++++++
    GroupAdvanced.add_argument('--max_cluster_bias_DEL',
                               help="Maximum distance to cluster read together for deletion.[%(default)s]",
                               default=1000,
                               type=int)
    GroupAdvanced.add_argument('--diff_ratio_merging_DEL',
                               help="Do not merge breakpoints with basepair identity more than [%(default)s] for deletion.",
                               default=0.5,
                               type=float)

    # ++++++INV++++++
    GroupAdvanced.add_argument('--max_cluster_bias_INV',
                               help="Maximum distance to cluster read together for inversion.[%(default)s]",
                               default=500,
                               type=int)

    # ++++++DUP++++++
    GroupAdvanced.add_argument('--max_cluster_bias_DUP',
                               help="Maximum distance to cluster read together for duplication.[%(default)s]",
                               default=500,
                               type=int)

    # ++++++TRA++++++
    GroupAdvanced.add_argument('--max_cluster_bias_TRA',
                               help="Maximum distance to cluster read together for translocation.[%(default)s]",
                               default=50,
                               type=int)
    GroupAdvanced.add_argument('--diff_ratio_filtering_TRA',
                               help="Filter breakpoints with basepair identity less than [%(default)s] for translocation.",
                               default=0.6,
                               type=float)

    args = parser.parse_args(argv)
    return args


def Generation_VCF_header(file, contiginfo, sample, argv):
    # General header
    file.write("##fileformat=VCFv4.2\n")
    file.write("##source=SKSV-%s/cuteSV-%s\n" % (VERSION_1,VERSION))
    import time
    file.write("##fileDate=%s\n" %
               (time.strftime('%Y-%m-%d %H:%M:%S %w-%Z', time.localtime())))
    for i in contiginfo:
        file.write("##contig=<ID=%s,length=%d>\n" % (i[0], i[1]))

    # Specific header
    # ALT
    file.write(
        "##ALT=<ID=INS,Description=\"Insertion of novel sequence relative to the reference\">\n")
    file.write(
        "##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">\n")
    file.write(
        "##ALT=<ID=DUP,Description=\"Region of elevated copy number relative to the reference\">\n")
    file.write("##ALT=<ID=INV,Description=\"Inversion of reference sequence\">\n")
    file.write("##ALT=<ID=BND,Description=\"Breakend of translocation\">\n")

    # INFO
    file.write(
        "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variant\">\n")
    file.write(
        "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">\n")
    file.write(
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    file.write(
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    file.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n")
    file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
    file.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n")
    file.write("##INFO=<ID=CILEN,Number=2,Type=Integer,Description=\"Confidence interval around inserted/deleted material between breakends\">\n")
    # file.write("##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n")
    file.write(
        "##INFO=<ID=RE,Number=1,Type=Integer,Description=\"Number of read support this record\">\n")
    file.write("##INFO=<ID=STRAND,Number=A,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">\n")
    file.write(
        "##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Supporting read names of SVs (comma separated)\">\n")
    file.write("##FILTER=<ID=q5,Description=\"Quality below 5\">\n")
    # FORMAT
    # file.write("\n")
    file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    file.write(
        "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# High-quality reference reads\">\n")
    file.write(
        "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# High-quality variant reads\">\n")
    file.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"# Phred-scaled genotype likelihoods rounded to the closest integer\">\n")
    file.write(
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"# Genotype quality\">\n")

    file.write("##CommandLine=\"cuteSV %s\"\n" % (" ".join(argv)))

