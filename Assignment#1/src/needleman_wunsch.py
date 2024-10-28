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

    for i in range(len_seq1 + 1):
        dp_matrix[i][0] = i * gap_penalty
    for j in range(len_seq2 + 1):
        dp_matrix[0][j] = j * gap_penalty

    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            pair = (seq1[i - 1], seq2[j - 1])
            score = substitution_matrix.get(pair, gap_penalty)
            match = dp_matrix[i - 1][j - 1] + score
            delete = dp_matrix[i - 1][j] + gap_penalty
            insert = dp_matrix[i][j - 1] + gap_penalty
            dp_matrix[i][j] = max(match, delete, insert)

    alignments = []
    def traceback(i, j, align1, align2):
        if len(alignments) >= n:
            return
        if i == 0 and j == 0:
            alignments.append((align1[::-1], align2[::-1]))
            return
        if i > 0 and dp_matrix[i][j] == dp_matrix[i-1][j] + gap_penalty:
            traceback(i-1, j, align1 + seq1[i-1], align2 + '-')
        if j > 0 and dp_matrix[i][j] == dp_matrix[i][j-1] + gap_penalty:
            traceback(i, j-1, align1 + '-', align2 + seq2[j-1])
        if i > 0 and j > 0:
            pair = (seq1[i-1], seq2[j-1])
            score = substitution_matrix.get(pair, gap_penalty)
            if dp_matrix[i][j] == dp_matrix[i-1][j-1] + score:
                traceback(i-1, j-1, align1 + seq1[i-1], align2 + seq2[j-1])