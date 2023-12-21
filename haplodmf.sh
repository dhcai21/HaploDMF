#! /bin/bash
### required arguments
file_sam=""
file_ref=""
### optional arguments
file_path='./result'
prefix="haplodmf"
thread=8
input_snv="None"
error_rate=0.1
signi_level=0.05
cond_pro=0.65
fre_snv=0.8
num_read_1=10
num_read_2=5
gap=15
smallest_snv=20
depth=5
mq=0
weight=0.2
lr=0.001
ep=20
alg=ward
c1=0.95
c2=0.90
li=20
batch_size=1024
s_pos=1
e_pos=10000000000
abundance=0


function help_info() {
	echo "Usage: $0 -i alignment.sam -r ref_genome.fasta [options]"
	echo ""
	echo "HaploDMF: viral Haplotypes reconstruction from long reads via Deep Matrix Factorization"
	echo ""
	echo "Author: Dehan CAI"
	echo "Date:   May 2022"
	echo ""
	echo "    -i | --input:                     alignment file (sam)"
	echo "    -r | --refernece:                 reference genome (fasta)"
	echo ""
	echo "    Options:"
	echo "    -o  | --out:                      Path where to output the results. (default:./result)"
	echo "    -p  | --prefix STR:               Prefix of output file. (default: rvhaplo)"
	echo "    -t  | --thread INT:               Number of CPU cores for multiprocessing. (default:8)"
	echo "    -is | --input_snv STR:            A file containing a set of SNV sites. (default:None)"
	echo "    -mq | --map_qual INT:             Smallest mapping quality in the alignment file. (default:0) "
	echo "    -sp | --start_pos INT:            Starting position for reconstructing haplotypes (default: 1)"
	echo "    -ep | --end_pos INT:              Ending position for reconstructing haplotypes (default: 1e10)"
	echo "    -a  | --abundance FLOAT:          Filter haplotypes with abundance less than a threshold (default: 0)"
	echo "    -h  | --help :                    Print help message."
	echo ""
	echo "        SNV Detection         "
	echo "    -e  | --error_rate FLOAT:         Sequencing error rate. (default: 0.1)"
	echo "    -s  | --signi_level FLOAT:        Significance level for binomial tests. (default: 0.05)"
	echo "    -cp | --cond_pro FLOAT:           A threshold of the maximum conditional probability for SNV sites. (default: 0.65)"
	echo "    -f  | --fre_snv FLOAT:            The most dominant base' frequency at a to-be-verified site should >= fre_snv. (default: 0.80)"
	echo "    -n1 | --num_read_1 INT:           Minimum # of reads for calculating the conditional probability given one conditional site. (default:10)"
	echo "    -n2 | --num_read_2 INT:           Minimum # of reads for calculating the conditional probability given more than one conditional sites. (default: 5)"
	echo "    -g  | --gap INT:                  Minimum length of gap between SNV sites for calculating the conditional probability. (default:15)"
	echo "    -ss | --smallest_snv INT:         Minimum # of SNV sites for haplotype construction. (default:20)"
	echo ""
	echo "      DMF model training"
	echo "    -w  | --weight FLOAT:             Weight (between 0 and 1) for loss 2. (default:0.2) "
	echo "    -lr | --learning_rate FLOAT:      Learning rate for DMF. (default:0.001) "
	echo "    -epo| --epoch INT:                Epoch for training. (default:20) "
	echo "    -bs | --batch_size INT:           Batch size for training. (default:1024) "
	echo ""
	echo "      CLustering  and reconstruction   "
	echo "    -al | --algorithm STR:            Algorithm for clustering: Hierarchical(ward) or Kmeans. (default:ward) "
	echo "    -d  | --depth INT:                Depth limitation for consensus sequences generated from clusters. (default:5) "
	echo "    -c1 | --cluster_thres1 STR:       Threshold 1 of edit distance with the increasing number of clusters. (default:0.95) "
	echo "    -c2 | --cluster_thres2 STR:       Threshold 2 of edit distance with the increasing number of clusters. (default:0.90) "
	echo "    -li | --largest_iteration INT:    Largest iteration for clustering. (default:20) "
	echo ""
	echo "    For further details of above arguments, please refer to https://github.com/dhcai21/HaploDMF   "
	echo ""
	exit 1
}

if [[ "$1" == "" ]];then
	help_info
	exit 1
fi

while [[ "$1" != "" ]]; do
	case "$1" in
		-h | --help ) ## print help message
		help_info
		exit 1
		;;
		-i | --input ) ### input sam file
		case "$2" in 
		"" )
			echo "Error: no sam file as input"
			exit 1
			;;
		*)
			file_sam="$2"
			if [[ "${file_sam:0:1}" == "-" ]]
			then
				echo "Error: no sam file as input"
				exit 1
			fi
			shift 2
			;;
		esac
		;;
		-r | --ref_genome) ### input reference genome
		case "$2" in 
		"")
			echo "Error: no fasta file as input"
			exit 1
			;;
		*)
			file_ref="$2"
			if [[ ""${file_ref:0:1}"" == "-" ]]
			then
				echo "Error: no fasta file as input"
				exit 1
			fi
			shift 2
			;;
		esac
		;;
		-o | --out )  ### output path
		case "$2" in 
		"" )
			echo "Error: no output path"
			exit 1
			;;
		*)
			file_path="$2"
			if [[ "${file_sam:0:1}" == "-" ]]
			then
				echo "Error: no output path"
				exit 1
			fi
			shift 2
			;;
		esac
		;;
		-p | --prefix )  ### prefix
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			prefix="$2"
			shift 2
			;;
		esac
		;;
		-t | --thread )  ### threads
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			thread="$2"
			shift 2
			;;
		esac
		;;
		-is | --input_snv )  ### input a file containing SNV sites from other tools
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			input_snv="$2"
			shift 2
			;;
		esac
		;;
		-mq | --map_qual )  ### depth limitation
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			mq="$2"
			shift 2
			;;
		esac
		;;
		-a | --abundance )  ### threshold of abundance
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			abundance="$2"
			shift 2
			;;
		esac
		;;
		-e | --error_rate )  ### error_rate
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			error_rate="$2"
			shift 2
			;;
		esac
		;;
		-s | --signi_level )  ### significance_level
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			signi_level="$2"
			shift 2
			;;
		esac
		;;
		-c | --cond_pro )  ### conditional_probability
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			cond_pro="$2"
			shift 2
			;;
		esac
		;;
		-f | --fre_snv )  ### determine the set of to-be-verified SNV sites
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			fre_snv="$2"
			shift 2
			;;
		esac
		;;
		-n1 | --num_read_1 )  ### number of reads for p(ai|aj)
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			num_read_1="$2"
			shift 2
			;;
		esac
		;;
		-n2 | --num_read_2 )  ### number of reads for p(ai|aj1,aj2,...,ajp)
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			num_read_2="$2"
			shift 2
			;;
		esac
		;;
		-g | --gap )  ### Minimum distance between SNV sites
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			gap="$2"
			shift 2
			;;
		esac
		;;
		-ss | --smallest_snv )  ### Minimum number of SNV sites for haplotype reconstruction
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			smallest_snv="$2"
			shift 2
			;;
		esac
		;;
		-w | --weight )  ### weight for loss
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			weight="$2"
			shift 2
			;;
		esac
		;;
		-lr | --learing_rate )  ### Learning rate
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			lr="$2"
			shift 2
			;;
		esac
		;;
		-epo | --epoch )  ### epoch
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			ep="$2"
			shift 2
			;;
		esac
		;;
		-bs | --batch_size )  ### batch size
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			batch_size="$2"
			shift 2
			;;
		esac
		;;
		-al | --algorithm )  ### clustering algorithm
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			alg="$2"
			shift 2
			;;
		esac
		;;
		-d | --depth )  ### depth limitation
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			depth="$2"
			shift 2
			;;
		esac
		;;
		-c1 | --cluster_thres1 )  ### Threshold 1 for clustering
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			c1="$2"
			shift 2
			;;
		esac
		;;
		-c2 | --cluster_thres2 )  ### Threshold 2 for clustering
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			c2="$2"
			shift 2
			;;
		esac
		;;
		-li | --largest_iteration )  ### largest iteration for clustering
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			li="$2"
			shift 2
			;;
		esac
		;;
		-sp | --start_pos )  ### Starting position
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			s_pos="$2"
			shift 2
			;;
		esac
		;;
		-ep | --end_pos )  ### Ending position
		case "$2" in 
		"" )
			echo "Error: no input for $1"
			exit 1
			;;
		*)
			e_pos="$2"
			shift 2
			;;
		esac
		;;
		*)
			echo "Error: unknow parameter $1"
			exit 1
	esac
done

if [[ "$file_sam" == "" ]];then
	echo "Error: no sam file input"
	exit 1
fi

if [[ "$file_ref" == "" ]];then
	echo "Error: no reference genome input"
	exit 1
fi

if [[ ${file_path:0-1:1} == "/" ]];then
	path_len=`expr ${#file_path}-1`
	file_prefix=$file_path$prefix
	file_path=${file_path:0:path_len}
else
	file_prefix=$file_path"/"$prefix
fi

##########  count nucleotide occurrence  ##########
echo "count nucleotide occurrence"
if [[ "$file_path" != "." ]];then
	rm -rf $file_path
	mkdir $file_path
fi
rm -rf $file_path"/alignment"
mkdir $file_path"/alignment"
file_len=`expr ${#file_sam}-4`
unique_sam=$file_path"/alignment/"$prefix".sam"
samtools view -h -F 0x900 -q $mq $file_sam > $unique_sam
file_bam=$file_path"/alignment/"$prefix".bam"
samtools view -b -S $unique_sam > $file_bam
rm $unique_sam
file_bam_sorted=$file_path"/alignment/"$prefix"_sorted.bam"
samtools sort $file_bam -o $file_bam_sorted
samtools index $file_bam_sorted
file_acgt=$file_prefix"_acgt.txt"
python ./src/count_frequency.py $file_bam_sorted $file_acgt
file_snv=$file_prefix"_snv.txt"

##########   SNV sites  ##########
if [[ "$input_snv" == "None" ]];then
	echo "SNV detection"
	python ./src/two_binomial.py $error_rate $signi_level $file_acgt $file_snv $thread $s_pos $e_pos
	## judge number of detected SNV sites
	size="$(wc -l <"$file_snv")"
	size="${size:0-1:1}"
	if (( $size != 0 ));then
		python ./src/out_haplotypes.py $prefix $file_bam_sorted $file_path 1 $s_pos $e_pos $abundance
		python ./src/extract_reads.py $file_path $prefix 1
		python ./src/run_medaka.py $file_path $prefix 1
		exit 1
	fi
else
	cp $input_snv $file_prefix"_snv.txt"
fi

##########   filter fake SNV sites and construct a frequency matrix   ##########

python ./src/fre_matrix.py $file_bam_sorted $file_snv $file_prefix $cond_pro $smallest_snv $num_read_1 $num_read_2 $gap $fre_snv $thread $input_snv

######### check the number of SNV sites ########
size="$(wc -l <"$file_snv")"
size="${size:0-1:1}"
if (( $size != 0 ));then
	python ./src/out_haplotypes.py $prefix $file_bam_sorted $file_path 1 $s_pos $e_pos $abundance
	python ./src/extract_reads.py $file_path $prefix 1
	python ./src/run_medaka.py $file_path $prefix 1
	exit 1
fi

##########   training and clustering  ##########
echo "Deep matrix factorization training and clustering"
python ./src/dmf.py --i $file_prefix"_matrix.pickle" --o $file_path --p $prefix --w $weight --l $lr --n $ep --algorithm $alg --c1 $c1 --c2 $c2 --largest_num $li --batch_size $batch_size

##########   reconstruct haplotypes  ##########
rm -rf $file_path"/clusters"
mkdir $file_path"/clusters"

echo "haplotypes reconstruction"

python ./src/out_haplotypes.py $prefix $file_bam_sorted $file_path x $s_pos $e_pos $abundance

echo "haplotypes polishment (medaka)"
python ./src/extract_reads.py $file_path $prefix x
python ./src/run_medaka.py $file_path $prefix x


rm $file_prefix"_matrix.pickle"
rm $file_prefix"_clusters.pickle"
rm -rf $file_path/medaka
echo "Finished"

exit 1
