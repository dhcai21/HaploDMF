# HaploDMF
viral Haplotyps reconstruction from long reads via Deep Matrix Factorization


### E-mail: dhcai2-c@my.cityu.edu.hk

### Dependencies:
* Conda
* Python >=3.8.13
* samtools >= 1.4.1
* [Medaka](https://github.com/nanoporetech/medaka)
* Pytorch>=1.10.0
* Required python package: pandas >= 0.18.0 tqdm, scipy

#### Install Dependencies
`conda create -n haplodmf python=3.8.13`<BR/>
`conda activate haplodmf`<BR/>
`conda install -c bioconda -c conda-forge medaka`<BR/>
`conda install pytorch=1.10.0=py3.8_cuda11.1_cudnn8.0.5_0 torchvision torchaudio cudatoolkit=11.1 -c pytorch`<BR/>
`conda install mcl`<BR/>
`pip install sklearn pandas tqdm`<BR/>
####
##### Possible problem
`'../lib/libcrypto.1.0.0.dylib' (no such file) when using samtools`

You can use the command:

`ln -s your_conda/rvhaplo/lib/libcrypto.your_exisiting_version.dylib your_conda/rvhaplo/lib/libcrypto.1.0.0.dylib`

## Usage
#### Initialization
`cd RVHaplo-main`<BR/>
`chmod +x rvhaplo.sh`
#### Command
`Example:   ./rvhaplo.sh -i alignment.sam -r reference.fasta -o result -p prefix -t 8`<BR/>

`If you want to output the SNV sites only:   ./rvhaplo.sh -i alignment.sam -r reference.fasta -os 1`<BR/>


`If you have multiple CPU cores, you can set the value of "-t" to accelerate the running.`<BR/>

```
required arguments:
    -i  | --input:                     alignment file (sam)
    -r  | --refernece:                 reference genome (fasta)

optional arguments:
    -h  | --help :                     Print help message.
    -o  | --out :                      Path where to output the results. (default: ./result)
    -p  | --prefix STR :               Prefix of output file. (default: rvhaplo)
    -t  | --thread INT :               Number of CPU cores for multiprocessing. (default: 8)
    -mq | --map_qual INT :             Smallest mapping quality for reads. (default: 0)
    -e  | --error_rate FLOAT :         Sequencing error rate. (default: 0.1)
    -s  | --signi_level FLOAT :        Significance level for binomial tests. (default: 0.05)
    -c  | --cond_pro FLOAT :           A threshold of the maximum conditional probability for SNV sites. (default: 0.65)
    -f  | --fre_snv :                  The most dominant base' frequency at a to-be-verified site should >= fre_snv. (default: 0.80)    
    -n1 | --num_read_1 INT :           Minimum # of reads for calculating the conditional probability given one conditional site. (default: 10)
    -n2 | --num_read_2 INT :           Minimum # of reads for calculating the conditional probability given more than one conditional sites. (default: 5)
    -g  | --gap INT :                  Minimum length of gap between SNV sites for calculating the conditional probability. (default: 15)
    -ss | --smallest_snv INT :         Minimum # of SNV sites for haplotype construction. (default: 20)
    -os | --only_snv (0 or 1) :        Only output the SNV sites without running the haplotype reconstruction part. (default: 0)
    -or | --overlap_read INT :         Minimum length of overlap for creating edges between two read in the read graph. (default: 5)
    -wr | --weight_read FLOAT :        Minimum weights of edges in the read graph. (default: 0.85)
    -m  | --mcl_inflaction FLOAT :     Inflaction of MCL algorithm. (default:2)
    -l  | --lar_cluster INT :          A threshold for seperating clusters into two groups based on sizes of clusters. (default: 50)
    -oc | --overlap_cluster INT :      A parameter related to the minimum overlap between consensus sequences of clusters. (default: 10)
    -d  | --depth INT :                Depth limitation for consensus sequences generated from clusters. (default: 5)
    -wc | --weight_cluster FLOAT :     Minimum weights between clusters in the hierarchical clustering (default: 0.8)
    -sp | --start_pos INT:             Starting position for generating consensus sequences (default: 1)
    -ep | --end_pos INT:               Ending position for generating consensus sequences. (default: 1e10)
    -a  | --abundance FLOAT :          A threshold for filtering low-abundance haplotypes. (default: 0.005)
```
`-t  | --thread`

If you have multiple CPU cores for running the tool, please use this parameter for accelerating the process. The default value is 8 (8 CPU cores).

`-mq  | --map_qual`

If your want to filter some reads with small mapping qualities, you can use this parameter (e.g., 20). When you data have a large number of reads, you can use this parameter to remove some bad-quality reads. This will help accelerate the running.

`-e  | --error_rate`

The sequencing error rate here can be roughly estimated. It will not significantly change the result for the bias between it and the ground truth. And we use 0.1 as the general sequencing error rate for TGS data.

`-s  | --signi_level`

Using a small significance level value may improve the precision of detected SNV sites obtained from binomial tests. But the small significance level value may reduce recall. Thus, we suggest using the default value 0.05.

`-c  | --cond_pro`

A threshold of the maximum conditional probability for SNV sites. If the maximum conditional probability of an SNV site is less than the threshold, this site will be recognized as a fake SNV site.

`-f  | --fre_snv`

Usually, sites containing fake SNVs caused by sequencing errors still have high frequencies of the most dominant bases. And those sites with small frequencies of the most dominant base are highly possible to contain real SNVs. Thus, we only verify part of potential sites obtained from the second binomial test to accelerate the verified process.

`-n1 | --num_read_1`

Minimum number of reads for calculating the conditional probability given one conditional site. For example, P(A|B).

`-n2 | --num_read_2`

The minimum number of reads for calculating the conditional probability given more than one conditional site. For example, P(A|B1,B2,B3,...). As the number of reads covering more SNVs sites will reduce, we allow a smaller number of reads for calculating the conditional probability given more conditional sites compared to only given one conditional site.

`-g  | --gap`

Because sites in close proximity tend to weaken the independence of sequencing errors, the distance (gap) between the target site and the given sites should be above a threshold.

`-ss | --smallest_snv`

As a small number of detected SNV sites indicates that only one strain in the sample is highly possible, the haplotype reconstruction process will stop and only output the detected SNV sites.

`-os | --only_snv`

Only output the SNV sites without running the haplotype reconstruction part.

`-or | --overlap_read`

Minimum length of overlap for creating edges between two read in the read graph.

`-wr | --weight_read`

Minimum weights of edges in the read graph.

`-m  | --mcl_inflaction`

The parameter "Inflaction" of the graph clustering algorithm Markov Cluster (MCL). Usually, using the default value 2 is enough here. For further details, please refer to https://micans.org/mcl/.

`-l  | --lar_cluster`

A threshold for seperating clusters into two groups based on sizes of clusters. Clusters with a larger size threshold will lead to more accurate consensus sequences but may miss some consensus sequences that are bridges between other consensus sequences. Users can modify the threshold based on the input number of reads. We suggest using one of (20, 30, 40, 50). 

`-oc | --overlap_cluster`

A parameter related to the minimum overlap between consensus sequences of clusters. The minimum overlap between consensus sequences is min(0.1 * x, oc), where x is the number of detected SNV sites.

`-d  | --depth`

Depth limitation for consensus sequences generated from clusters. The total number of bases (A, C, G, T) at a site should be larger than the threshold. Otherwise, ignore this site (use '-').

`-wc | --weight_cluster`

Minimum weights between clusters in the hierarchical clustering.

`-sp  | --start_pos`

Starting position for generating consensus sequences. (1-index)

`-ep  | --end_pos`

Ending position for generating consensus sequences. A large default value is for covering the whole genome. (1-index)

A threshold for filtering low-abundance haplotypes.

`-a  | --abundance`

A threshold for filtering low-abundance haplotypes.



## Output Results
All reconstructed haplotypes are summarized in a file "*_haplotypes.fasta". Below is an example of three haplotypes.
```
>haplotype_0_length_9181_abundance_0.50
GGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGGTCTCTGGCTAACTAGGGAACC...
>haplotype_1_length_9178_abundance_0.30
GTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCC...
>haplotype_2_length_9180_abundance_0.20
GGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGGACC...
```

#### Clusters of Reads
For the corresponding read clusters of haplotypes, the alignment file of each cluster can be found in the folder "clusters".

#### SNV sites
You can find out the SNV sites (1-index on the reference genome) in the file "*_SNV.txt".

The base (A, C, G, T, etc.) counts from the alignments between reads and the reference genome can be found in the file "*_acgt.txt".

`
If you want to output the SNV sites only without running the haplotype reconstruction part, use the command "-os 1" or "--only_snv 1".
`

## Simulation data
Simulation datasets can be downloaded from https://drive.google.com/drive/folders/16azUqV6thGJyBThR0OaGfNgjXYhD_lEE?usp=sharing

## Tested data
The tested data can be downloaded from https://drive.google.com/drive/folders/1A4I722MO4Ujnvu_B3pugq6v6N4t1koHd?usp=sharing

`./rvhaplo.sh -i test.sam -r reference.fasta`<BR/>


