# -*- coding: utf-8 -*-
"""transcription-binding-predictor.ipynb
"""

'''
goal is to predict whether a sequence of DNA will be bound by a specific transcription 
factor in a given condition given a PWM (postional weight matrix)
'''

import argparse
import numpy as np

# Create the argument parser
parser = argparse.ArgumentParser(description='File paths')

# Add the file path arguments
parser.add_argument('test_file', type=str, help='Path to the test file')
parser.add_argument('PWM_file', type=str, help='Path to the PWM file')

# Parse the command-line arguments
args = parser.parse_args()

#function to read in fasta files
def read_fasta_file(filename):
    with open(filename, 'r') as f:
        read_id = None
        read_seq = ""
        for line in f:
            if line.startswith(">"):
                if read_id is not None:
                    yield (read_id, read_seq)
                read_id = line.strip()[1:]
                read_seq = ""
            else:
                read_seq += line.strip()
        if read_id is not None:
            yield (read_id, read_seq)

#function to read in PWM
def read_profile_matrix(file_path):
    profile_matrix = np.loadtxt(file_path, delimiter='\t')
    return profile_matrix

#function to find the reverse complement of a sequence
def ReverseComplement(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement_seq = [complement_dict[base] for base in reversed(seq)]
    complement = ''.join(complement_seq)
    return complement

#function to calculate proabilities for kmers 
def calculate_probability(kmer, pwm):
    seq1 = 'ACGT0123'
    seq_dict = {seq1[i]: i for i in range(4)}
    p = 1
    for i in range(len(kmer)):
        nucleotide = kmer[i]
        index = seq_dict[nucleotide]
        p *= pwm[index][i]
    return p

#function to keep the highest probability among a kmer in a sequence 
def find_max_probability(sequence, kmer_length, pwm):
    max_prob = 0.0

    for i in range(len(sequence) - kmer_length + 1):
        kmer = sequence[i:i+kmer_length]
        probability = calculate_probability(kmer, pwm)
        if probability > max_prob:
            max_prob = probability

    return max_prob

#function to compute the highest probability for each forward sequence 
def all_forward_probabilites(test_reads, kmer_length, pwm):
    all_forward_probs = []
    
    for read_id, read_seq in test_reads:
        max_prob = find_max_probability(read_seq, kmer_length, pwm)
        read_probs = {'read_id': read_id, 'max_prob': max_prob}
        all_forward_probs.append(read_probs)
    
    return all_forward_probs


#function to compute the highest probability for each reverse sequence 
def all_reverse_probabilites(test_reads, kmer_length, pwm):
    all_rev_probs = []
    
    for read_id, read_seq in test_reads:
        rev_read = ReverseComplement(read_seq)
        max_prob = find_max_probability(rev_read, kmer_length, pwm)
        read_probs = {'read_id': read_id, 'max_prob': max_prob}
        all_rev_probs.append(read_probs)
    
    return all_rev_probs


# Function to merge and write the top 2000 sequences with the highest probabilities to a file
def merge_and_write_top_sequences(forward_probs, reverse_probs, output_file, top_k=2000):
    # Merge forward and reverse probabilities
    merged_probs = forward_probs + reverse_probs

    # Sort the sequences by max_prob in descending order
    merged_probs.sort(key=lambda x: x['max_prob'], reverse=True)

    # Write the top sequences to the output file
    written_seqs = set()  # Track sequences that have been written
    num_written_seqs = 0  # Track the number of written sequences

    with open(output_file, 'w') as file:
        i = 0
        while num_written_seqs < top_k and i < len(merged_probs):
            sequence = merged_probs[i]['read_id']
            if sequence not in written_seqs:
                file.write(f"{sequence}\n")
                written_seqs.add(sequence)
                num_written_seqs += 1
            i += 1


#call functions

# Access the file paths using argparse
test_file = args.test_file
PWM_file = args.PWM_file

test_reads = read_fasta_file(test_file)
PWM = read_profile_matrix(PWM_file)
kmer_length = 21

output_file = 'predictions.txt'

for_probs = all_forward_probabilites(test_reads, kmer_length, PWM) 

test_reads = read_fasta_file(test_file)
rev_probs = all_reverse_probabilites(test_reads, kmer_length, PWM) 

# Output a predictions.txt file of the top 2000 sequences that have the highest probability of binding to a specific transcription factor
merge_and_write_top_sequences(for_probs, rev_probs, output_file, top_k=2000)
