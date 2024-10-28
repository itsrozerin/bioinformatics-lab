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
    
