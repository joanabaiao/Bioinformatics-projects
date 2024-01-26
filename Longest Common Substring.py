# LONGEST COMMON SUBSTRING

import numpy as np

seq1 = "ACGTCATCA"
seq2 = "TAGTGTCA"  

m = len(seq1)
n = len(seq2)

matrix = np.zeros((m+1, n+1))
longest = 0
LCS = ""

for i in range(m):
    for j in range(n):

        if seq1[i] == seq2[j]:
            c = matrix[i][j] + 1
            matrix[i+1][j+1] = c

            if c >= longest:
                longest = c
                LCS = seq1[int(i-c+1):i+1]

print(matrix)
print()
print("Longest common substring:", LCS)
print("Length of the substring:", int(longest))

