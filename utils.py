import argparse
import os
import sys
import numpy as np


def cmd():
    parser = argparse.ArgumentParser(description="ARGUMENTS")
    
    # Input data
    parser.add_argument(
        "--i",
        default='frequency_matrix.pickle',
        type=str,
        help="frequency_matrix (default: frequency_matrix.pickle)")
    
    # Output data
    parser.add_argument(
        "--o",
        default='result',
        type=str,
        help="output dir (default: result)")
    
    # Prefix of output file
    parser.add_argument(
        "--p",
        default='haplodmf',
        type=str,
        help="file prefix (default: haplodmf)")
    
    # Weight of losses
    parser.add_argument(
        "--w",
        default=0.2,
        type=float,
        help="weight of two losses (default: 0.2)")
    
    # for network training
    parser.add_argument(
        '--l',
        default=0.001,
        type=float,
        help='Initial learning rate (default: 0.001)'
    )
    
    # for network training
    parser.add_argument(
        '--n',
        default=20,
        type=int,
        help="Number of epochs to train (default: 20)"
    )
    
    # for clustering
    parser.add_argument(
        '--algorithm',
        type=str,
        default="ward",
        help="Custering algorithm: ward or kmeans. (default: ward)"
    )
    
    parser.add_argument(
        '--c1',
        type=float,
        default=0.95,
        help="first threshold for the change of edit distance with the increase of haplotypes' number. (default: 0.95)"
    )

    parser.add_argument(
        '--c2',
        type=float,
        default=0.90,
        help="second threshold for the change of edit distance with the increase of haplotypes' number. (default: 0.90)"
    )
    
    parser.add_argument(
        '--largest_num',
        type=int,
        default=20,
        help="Largest number of haplotypes. (default: 20)"
    )
    
    parser.add_argument(
        '--batch_size',
        type=int,
        default=1024,
        help="Batch size. (default: 1024)"
    )

    args = parser.parse_args()

    return args
