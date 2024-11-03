import numpy as np

# Traceback function to extract optimal alignments
def trace_back(dp, seq1, seq2, substitution_matrix, gap_penalty, max_alignments, start_pos=None):
    alignments = []
    m, n = len(seq1), len(seq2)
    
    # Set the starting position for Smith-Waterman (highest score cell)
    if start_pos:
        i, j = start_pos
        score = dp[i][j]  # Start with the highest score
    else:
        i, j = m, n  # Start from the bottom-right corner for Needleman-Wunsch
        score = dp[m][n]  # Take the bottom-right score for global alignment

    # Recursive function for traceback
    def traceback_recursive(i, j, aligned_seq1, aligned_seq2):
        if len(alignments) >= max_alignments:
            return
        if dp[i][j] == 0:  # Stop at zero for Smith-Waterman
            alignments.append((aligned_seq1[::-1], aligned_seq2[::-1], score))
            return
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + substitution_matrix[(seq1[i-1], seq2[j-1])]:
            traceback_recursive(i-1, j-1, aligned_seq1 + seq1[i-1], aligned_seq2 + seq2[j-1])
        if i > 0 and dp[i][j] == dp[i-1][j] + gap_penalty:
            traceback_recursive(i-1, j, aligned_seq1 + seq1[i-1], aligned_seq2 + '-')
        if j > 0 and dp[i][j] == dp[i][j-1] + gap_penalty:
            traceback_recursive(i, j-1, aligned_seq1 + '-', aligned_seq2 + seq2[j-1])
    
    traceback_recursive(i, j, "", "")
    print("Traceback results:", alignments)
    return alignments