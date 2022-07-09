# HaploDMF
viral Haplotyps reconstruction from long reads via Deep Matrix Factorization


### E-mail: dhcai2-c@my.cityu.edu.hk

### Dependencies:
* Conda
* Python >=3.8.13
* samtools
* [Medaka](https://github.com/nanoporetech/medaka)
* Pytorch>=1.10.0
* Required python package: pysam, sklearn, pandas, tqdm, scipy

### An easiler way to install
After cloning this respository, you can use anaconda to install the **haplodmf.yaml**. This will install all packages you need with gpu mode (make sure you have installed cuda on your system to use the gpu version. Othervise, it will run with cpu version). The command is: `conda env create -f haplodmf.yaml -n haplodmf`

* For cpu version pytorch: `conda install pytorch torchvision torchaudio cpuonly -c pytorch`
* For gpu version pytorch: Search [pytorch](https://pytorch.org/) to find the correct cuda version according to your computer

#### An optional way to install
`conda create -n haplodmf python=3.8.13`<BR/>
`conda activate haplodmf`<BR/>
`conda install -c bioconda -c conda-forge medaka`<BR/>
`conda install pytorch=1.10.0=py3.8_cuda11.1_cudnn8.0.5_0 torchvision torchaudio cudatoolkit=11.1 -c pytorch`<BR/>
`pip install sklearn pandas tqdm`<BR/>


##### Possible problems
* '../lib/libcrypto.1.0.0.dylib' (no such file) when using samtools

  You can use the command:

  `ln -s your_conda/rvhaplo/lib/libcrypto.your_exisiting_version.dylib your_conda/rvhaplo/lib/libcrypto.1.0.0.dylib`

* AttributeError: module 'distutils' has no attribute 'version'

  You can use the command:

  `pip install setuptools==59.5.0`

## Usage
#### Initialization
`cd HaploDMF-main`<BR/>
`chmod +x haplodmf.sh`
#### Command
`Example:   ./haplodmf.sh -i alignment.sam -r reference.fasta -o result -p prefix`<BR/>

```
Required arguments:
    -i | --input:                     alignment file (sam)
    -r | --refernece:                 reference genome (fasta)

Options:
    -o  | --out:                      Path where to output the results. (default:./result)
    -p  | --prefix STR:               Prefix of output file. (default: rvhaplo)
    -is | --input_snv STR:            A file containing a set of SNV sites. (default:None)
    -mq | --map_qual INT:             Smallest mapping quality in the alignment file. (default:0) 
    -sp | --start_pos INT:            Starting position for reconstructing haplotypes (default: 1)
    -ep | --end_pos INT:              Ending position for reconstructing haplotypes (default: 1e10)
    -h  | --help :                    Print help message.

    SNV Detection:
    -t  | --thread INT:               Number of CPU cores for multiprocessing. (default:8)
    -e  | --error_rate FLOAT:         Sequencing error rate. (default: 0.1)
    -s  | --signi_level FLOAT:        Significance level for binomial tests. (default: 0.05)
    -cp | --cond_pro FLOAT:           A threshold of the maximum conditional probability for SNV sites. (default: 0.65)
    -f  | --fre_snv FLOAT:            The most dominant base' frequency at a to-be-verified site should >= fre_snv. (default: 0.80)
    -n1 | --num_read_1 INT:           Minimum # of reads for calculating the conditional probability given one conditional site. (default:10)
    -n2 | --num_read_2 INT:           Minimum # of reads for calculating the conditional probability given more than one conditional sites. (default: 5)
    -g  | --gap INT:                  Minimum length of gap between SNV sites for calculating the conditional probability. (default:15)
    -ss | --smallest_snv INT:         Minimum # of SNV sites for haplotype construction. (default:20)

    DMF model training
    -w  | --weight FLOAT:             Weight (between 0 and 1) for loss 2. (default:0.2) 
    -lr | --learning_rate FLOAT:      Learning rate for DMF. (default:0.001) 
    -epo| --epoch INT:                Epoch for training. (default:20) 
    -bs | --batch_size INT:           Batch size for training. (default:1024) 

    Clustering  and reconstruction   
    -al | --algorithm STR:            Algorithm for clustering: Hierarchical(ward) or Kmeans. (default:ward) 
    -d  | --depth INT:                Depth limitation for consensus sequences generated from clusters. (default:5) 
    -c1 | --cluster_thres1 STR:       Threshold 1 of edit distance with the increasing number of clusters. (default:0.95) 
    -c2 | --cluster_thres2 STR:       Threshold 2 of edit distance with the increasing number of clusters. (default:0.90) 
    -li | --largest_iteration INT:    Largest iteration for clustering. (default:20) 
```

`-mq  | --map_qual`

If you want to filter some reads with small mapping qualities, you can use this parameter (e.g., 20). When you data have a large number of reads, you can use this parameter to remove some bad-quality reads. This will help accelerate the running.

`-is  | --input_snv`

You can use other SNV detection tools to obtain a set of SNV sites (1-index on a reference genome). And then input a file containing the SNV sites to HaploDMF for downstream analysis. An example of the file can be found in the folder *tested_data*.


`-sp  | --start_pos`

Starting position for reconstructing haplotypes. (1-index)

`-ep  | --end_pos`

Ending position for reconstructing haplotypes. A large default value is for covering the whole genome. (1-index)


## Run HaploDMF on tested data
`./haplodmf.sh -i test.sam -r reference.fasta -o test_result -p test`<BR/>

## Output Results
The trained model can be found in the folder "log".

All reconstructed haplotypes (Polished by Medaka) are summarized in a file "*_haplotypes.fasta". Below is an example of three haplotypes.

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


## Simulation data
Simulation datasets can be downloaded from 


