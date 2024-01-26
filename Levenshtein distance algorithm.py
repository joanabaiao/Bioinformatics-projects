# LEVENSHTEIN DISTANCE ALGORITHM

import numpy as np

seq1 = "sunday"
seq2 = "saturday"

m = len(seq1) + 1
n = len(seq2) + 1

matrix = np.zeros([m, n])

for i in range(m):
    matrix[i, 0] = i
for j in range(n):
    matrix[0, j] = j

for i in range(1, m):
    for j in range(1, n):

        cost = 0
        if seq1[i - 1] != seq2[j - 1]:
            cost = 1
        
        score_up = matrix[i-1][j] + 1  # deletion
        score_left = matrix[i][j-1] + 1  # insertion
        score_diagonal = matrix[i-1][j-1] + cost  # mutation or conservation

        
        matrix[i][j] = min(score_up, score_left, score_diagonal)

print(matrix, "\n")
print("Minimum number:", matrix[m-1][n-1])
