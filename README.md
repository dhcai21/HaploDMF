# HaploDMF
HaploDMF: viral Haplotype reconstruction from long reads via Deep Matrix Factorization


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
```
conda create -n haplodmf python=3.8.13

conda activate haplodmf

conda install -c bioconda -c conda-forge medaka

conda install pytorch torchvision torchaudio cudatoolkit -c pytorch  (Install a correct version of pytorch according to your computer)

pip install sklearn pandas tqdm
```

##### Possible problems
* '../lib/libcrypto.1.0.0.dylib' (no such file) when using samtools

  You can use the command:

  `ln -s your_conda/rvhaplo/lib/libcrypto.your_exisiting_version.dylib your_conda/rvhaplo/lib/libcrypto.1.0.0.dylib`

* AttributeError: module 'distutils' has no attribute 'version'

  You can use the command:

  `pip install setuptools==59.5.0`

## Usage
#### Initialization
`cd HaploDMF`<BR/>
`chmod +x haplodmf.sh`
#### Command
`Example:   ./haplodmf.sh -i alignment.sam -r reference.fasta -o result -p prefix`<BR/>

#### Parameters

In HaploDMF, there are several groups of parameters with different purposes: parameters for detecting SNV sites, training the neural network, and clustering. Most of these parameters are provided for other developers who might want to modify our codes for other purposes.

Only the following parameters may affect the result:
`1. In the training process, we used loss2 to help learn latent features between distant reads. When the sequencing coverage is large, loss2 will not affect the result significantly. Thus, the parameter of controlling the weights between two losses is not sensitive when the sequencing coverage is high. When the coverage is low, loss2’s weight can be set to 0.1-0.2. `

`2. We used the “elbow method” to determine the number of clusters (haplotypes). Empirically, using the default values is suitable for most cases. Stringent values may output a few repetitive haplotypes with low abundance while relaxed values may miss the low-abundant haplotypes. But users can adjust the parameters according to the abundance of haplotypes, number of reads supporting haplotypes, etc., which are provided in the output file.`


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
    -a  | --abundance FLOAT:          Filter haplotypes with abundance less than a threshold (default: 0)
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
    -al | --algorithm STR:            Algorithm for clustering: Hierarchical(ward) or KMeans. (default:ward) 
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

`-a  | --abundance `

This parameter is to filter low-abundance haplotypes from the output file. If you use "-a 0.05", then the output file will only contain haplotypes with abundance >= 0.05 (i.e., 5%). The default setting will keep all the haplotypes in the file. And the output file contains information including the number of reads, depth, and abundance for each haplotype. Theoretically, the threshold of abundance is related to the error rate. However, it is hard to define a clear relationship between error rate and abundance. According to our experience, when the sequencing error rate is high (e.g., 18%), the confidence of haplotypes with an abundance of less than 1% is small. Users need to consider other information like the number of reads to determine whether the low-abundant haplotypes are reliable. 

`SNV detection`

Please refer to [RVHaplo](https://github.com/dhcai21/RVHaplo)

`-w  | --weight`

There are two losses in HaploDMF. Minimizing loss1 learns latent features of reads so that reads sharing more SNVs have similar latent features. Loss2 is to help learn latent features for distant SNV sites. If the distant SNV sites contain SNVs from some haplotypes, latent features of distant reads (without overlap) can be better generated. Because loss1 is the main loss to learn the latent features, we assign a large weight to it and a small weight to loss2. Here, -w/--weight is to assign the weight to loss2. According to our experiments, we suggest selecting a value between (0,0.2]. If the coverage of your data is large and without drop-off in the middle, this parameter will not affect the result significantly. 

`-lr  | --learning_rate`

The learning rate when training the neural network. 

`-epo  | --epoch`

The epoch for training. As the training process can converge fast, we use a small epoch here to reduce the training time. 

`-bs | --batch_size`

The batch size for training. If the memory of your device is small, you can reduce the batch size.

`-al | --algorithm`

The default algorithm for cluster reads is hierarchical clustering (Ward's method). Because hierarchical clustering is time-consuming when the dataset is large, we provide an alternative clustering algorithm "KMEAMS". According to our experiments, the clustering results between these two methods are similar.

`-d | --depth`

The depth limitation for generating consensus sequences. Because the consensus sequence contains more errors in the ending regions due to the small depth, we set a threshold of depth to remove these regions when calculating the distance between reads and consensus sequences for determining the number of clusters. As the depth is small, removing these regions will not affect the final number of clusters. 

`-c1 | --cluster_thres1`

We use the "elbow method" to determine the number of clusters. If the number of different bases between reads and their consensus sequences does not change significantly (e.g., > 0.95) with the increasing number of clusters, the iteration will stop. 

`-c2 | --cluster_thres2`

Because the number of different bases may not decrease significantly at the beginning iteration when the sequencing error rate is high, we use this threshold to double-check the trend of different bases with the increasing number of clusters. For further details, please refer to the supplementary file.

`-li | --largest_iteration`

Because the high values of "-c1" and "-c2" may not be able to stop the iteration, we set the largest iteration to avoid the infinite iteration. This value should be large than the number of haplotypes. Usually, the default value (20) is suitable for most datasets.

## Run HaploDMF on tested data
`./haplodmf.sh -i test.sam -r reference.fasta -o test_result -p test`<BR/>

## Output Results
The trained model can be found in the folder "log".

All reconstructed haplotypes (Polished by Medaka) are summarized in a file "*_haplotypes.fasta". Below is an example of three haplotypes.

```
>haplotype_0_length_9181_abundance_0.50_number_of_reads_5000_depth_500
GGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGGTCTCTGGCTAACTAGGGAACC...
>haplotype_1_length_9178_abundance_0.30_number_of_reads_3000_depth_300
GTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCC...
>haplotype_2_length_9180_abundance_0.20_number_of_reads_2000_depth_200
GGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGGACC...
```

#### Clusters of Reads
For the corresponding read clusters of haplotypes, the alignment file of each cluster can be found in the folder "clusters".

#### SNV sites
You can find out the SNV sites (1-index on the reference genome) in the file "*_SNV.txt".

The base (A, C, G, T, etc.) counts from the alignments between reads and the reference genome can be found in the file "*_acgt.txt".


## Simulation data
Simulation datasets can be downloaded at https://portland-my.sharepoint.com/:f:/g/personal/dhcai2-c_my_cityu_edu_hk/Eiey3SkmQXxGlEfWhpgzCVcBmxQ12kbaE6m-djT6c5rWRA?e=Gkz1Di 


