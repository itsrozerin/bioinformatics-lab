import numpy as np
from traceback import trace_back
from file_operations import write_alignments_to_file, read_substitution_matrix

# Needleman-Wunsch (NW) global alignment function
def needleman_wunsch(substitution_matrix, gap_penalty, max_alignments, nw_output_filename):
    seq1 = "TATA"  # First example DNA sequence
    seq2 = "ATAT"  # econd example DNA sequence
    substitution_matrix = read_substitution_matrix("dna_matrix.csv")
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m+1, n+1), dtype=int) # Here we need to intialize a matrix containing (m+1)*(n+1). 
    # Because we are going to have gap penalties in the first row and column.
    
    # Initializing the dynamic programming matrix with gap penalties. We need gap penalties in the first row and column.
    for i in range(1, m+1):
        dp[i][0] = dp[i-1][0] + gap_penalty
    for j in range(1, n+1):
        dp[0][j] = dp[0][j-1] + gap_penalty

    # Fill in the dynamic programming matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = dp[i-1][j-1] + substitution_matrix[(seq1[i-1], seq2[j-1])]
            delete = dp[i-1][j] + gap_penalty
            insert = dp[i][j-1] + gap_penalty
            dp[i][j] = max(match, delete, insert)
    
    # Find the optimal alignments and format the output
    alignments = trace_back(dp, seq1, seq2, substitution_matrix, gap_penalty, max_alignments)

    #To check if the dynamic programming matrix is created properly
    print(dp)
    write_alignments_to_file(alignments, nw_output_filename)

# Smith-Waterman (SW) local alignment function
def smith_waterman(substitution_matrix, gap_penalty, max_alignments, sw_output_filename):
    seq1 = "TATA"  # first example DNA sequence
    seq2 = "ATAT"  # second example DNA sequence
    substitution_matrix = read_substitution_matrix("dna_matrix.csv")
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m + 1, n + 1), dtype=int) # Here we need to intialize a matrix containing (m+1)*(n+1).
    # Beacuse we are going to have zeroes in the first column and row.
    
    max_score = 0 # Initialize max score as 0
    max_pos = (0, 0)
    
    # Fill in the DP matrix with local alignment rules
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i-1][j-1] + substitution_matrix[(seq1[i-1], seq2[j-1])]
            delete = dp[i-1][j] + gap_penalty
            insert = dp[i][j-1] + gap_penalty
            dp[i][j] = max(0, match, delete, insert)  # Zero for local alignment. We don't have values below zero.
            
            # Track the cell with the highest score
            if dp[i][j] > max_score:
                max_score = dp[i][j]
                max_pos = (i, j)
    
    # Perform traceback starting from the cell with the highest score
    alignments = trace_back(dp, seq1, seq2, substitution_matrix, gap_penalty, max_alignments, start_pos=max_pos)
    print(dp)
    write_alignments_to_file(alignments, sw_output_filename)
