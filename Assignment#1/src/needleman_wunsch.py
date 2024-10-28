import csv
import numpy as np

# Function to read the substitution matrix from a CSV file
def read_substitution_matrix(filepath):
    substitution_matrix = {} # Initialize an empty dictionary for the substitution matrix
    with open(filepath, mode = 'r') as file:
        reader = csv.reader(file) # Use csv reader to read the file row by row
        headers = next(reader)
        for i, row in enumerate(reader):
            pair = (headers[i+1], headers[j+1])
            substitution_matrix[pair] = int(score) # Convert score to integer and store
        return substitution_matrix
    

def needleman_wunsch(seq1, seq2, n, substitution_matrix_path, gap_penalty, output_file):
    # Read substitution matrix from provided file path
    substitution_matrix = read_substitution_matrix(substitution_matrix_path)
    len_seq1, len_seq2 = len(seq1), len(seq2)
    # Initialize dynamic programming (DP) matrix with zeros
    dp_matrix = np.zeros((len_seq1 + 1, len_seq2 + 1))

    # Initialize first row and column with gap penalties
    for i in range(len_seq1 + 1):
        dp_matrix[i][0] = i * gap_penalty
    for j in range(len_seq2 + 1):
        dp_matrix[0][j] = j * gap_penalty

    # Fill the DP matrix based on substitution scores and gap penalties
    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            pair = (seq1[i - 1], seq2[j - 1]) # Create pair for the current nucleotide comparison
            score = substitution_matrix.get(pair, gap_penalty) # Get score from matrix or default to gap penalty
            match = dp_matrix[i - 1][j - 1] + score
            delete = dp_matrix[i - 1][j] + gap_penalty
            insert = dp_matrix[i][j - 1] + gap_penalty
            dp_matrix[i][j] = max(match, delete, insert) # Use maximum score for optimal alignment

    
    # Initialize list to store multiple alignments
    alignments = []
     # Recursive function to trace back and build optimal alignments
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
            
        # Start traceback from bottom-right cell
    traceback(len_seq1, len_seq2, '', '')

   

    # Write results to output file
    with open(output_file, 'w') as f:
        for idx, (align1, align2) in enumerate(alignments[:n], 1):
            f.write(f"Global alignment no. {idx}:\n")
            f.write(f"{align1}\n{align2}\n")
            f.write(f"Score: {int(dp_matrix[len_seq1][len_seq2])}\n\n")