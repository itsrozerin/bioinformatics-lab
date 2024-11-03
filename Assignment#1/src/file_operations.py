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
