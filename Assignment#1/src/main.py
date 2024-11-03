"""
Finding out about similarities in different DNA sequences is an important task in Bioinformatics. 
In this project, I implemented two different DNA sequence allignment algortihms: Needleman-Wunsch and Smith-Waterman.
Learnt more about what can I do with Python and its many useful libraries like NumPy.
Also, shoutout to Needleman, Wunsch, Smith and Waterson for coming up with these brilliant ways of analizing sequences :)
Rozerin Dag, 2024
"""

import numpy 
from alignments import needleman_wunsch, smith_waterman

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