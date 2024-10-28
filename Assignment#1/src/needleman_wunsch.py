import csv
import numpy

def read_substitution_matrix(filepath):
    substitution_matrix = {}
    with open(filepath, mode = 'r') as file:
        reader = csv.reader(file)
        headers = next(reader)
        for i, row in enumerate(reader):
            pair = (headers[i+1], headers[j+1])
            substitution_matrix[pair] = int(score)
        return substitution_matrix
    

def needleman_wunsch(seq1, seq2, n, substitution_matrix_path, gap_penalty, output_file):
    substitution_matrix = read_substitution_matrix(substitution_matrix_path)
    len_seq1, len_seq2 = len(seq1), len(seq2)
    dp_matrix = np.zeros((len_seq1 + 1, len_seq2 + 1))

    