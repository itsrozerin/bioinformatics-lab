import numpy as np

# Function to read the substitution matrix from file
def read_substitution_matrix(filepath):
    matrix = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()
        headers = lines[0].strip().split(",")[1:]  # Separating the first row and first column as headers
        
        for i, line in enumerate(lines[1:]):  # We need to skip the header row
            values = line.strip().split(",")
            row = values[0]  # The first cell contains the row header
            for j, val in enumerate(values[1:]):  # Process all other cells except the first because it is the header cell
                col = headers[j]
                matrix[(row, col)] = int(val)  # Store values as integers

    print(matrix)   # To see if the matrix is created properly     
    return matrix

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

# Function to write alignments to a file
def write_alignments_to_file(alignments, output):
    if not alignments:
        print("No alignments found.")
        return
    with open(output, 'w') as f:
        for i, (alignment_seq1, alignment_seq2, score) in enumerate(alignments, 1):
            f.write(f"Alignment {i}:\n")
            f.write(f"{alignment_seq1}\n")
            f.write(f"{alignment_seq2}\n")
            f.write(f"Score: {score}\n\n")

# Main function
def main():

    matrix_filepath = "dna_matrix.csv"
    gap_penalty = -2
    max_alignments = 2
    nw_output_filename = "nw_alignments.txt"
    sw_output_filename = "sw_alignments.txt"
    
    # Run Needleman-Wunsch and Smith-Waterman algorithms
    nw_alignments = needleman_wunsch(matrix_filepath, gap_penalty, max_alignments,nw_output_filename)
    sw_alignments = smith_waterman(matrix_filepath, gap_penalty, max_alignments, sw_output_filename)    

if __name__ == "__main__":
    main()