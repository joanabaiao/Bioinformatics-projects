# LONGEST COMMON SUBSEQUENCE

import numpy as np

seq1 = "ACGTCATCA"
seq2 = "TAGTGTCA"
m = len(seq1) 
n = len(seq2) 

matrix = np.zeros((m+1, n+1))

for i in range(m+1): 
    for j in range(n+1): 
        if i == 0 or j == 0: 
            matrix[i][j] = 0
        elif seq1[i-1] == seq2[j-1]: 
            matrix[i][j] = matrix[i-1][j-1] + 1
        else: 
            matrix[i][j] = max(matrix[i-1][j], matrix[i][j-1]) 

print(matrix)
index = int(matrix[m][n])

string = [""] * (index+1) 
string[index] = "" 

i = m 
j = n 
while i > 0 and j > 0: 

    if seq1[i-1] == seq2[j-1]: 
        string[index-1] = seq1[i-1] 
        i = i - 1
        j = j - 1
        index = index - 1

    elif matrix[i-1][j] > matrix[i][j-1]: 
        i = i - 1
    else: 
        j = j - 1

LCSS = "".join(string)
print("\nLongest common subsequence:", LCSS)

###########################################################################################

print("\n-------------------------\n")

matrix = np.zeros((len(seq1)+1, len(seq2)+1))

for i in range(len(seq1)+1):
    matrix [i, 0] = 0
for j in range(len(seq2)+1):
    matrix [0, j] = 0
    
for i in range(1,len(seq1)+1):
    for j in range(1,len(seq2)+1):   

        if seq1[i-1] == seq2[j-1]:
            matrix[i][j] = matrix[i-1][j-1] + 1           
        else:
            matrix[i][j] = max(matrix[i-1,j], matrix[i,j-1], matrix[i-1,j-1])   

print(matrix)

def traceback(matrix):
    i = matrix.shape[0]-1
    j = matrix.shape[1]-1

    aligned_seq = []

    while i != 0 and j != 0:
        up = matrix[i-1, j]
        side = matrix[i, j-1]

        if seq1[i-1] == seq2[j-1]:
            aligned_seq.append(seq1[i-1])
            i = i-1
            j = j-1

        else:
            cell = max(up, side)

            if cell == up:
                #aligned_seq.append(" ")
                i = i-1

            elif cell == side:
                #aligned_seq.append(" ")
                j = j-1

    aligned_seq.reverse()
    aligned_seq = "".join(aligned_seq)

    return aligned_seq

aligned_seq = traceback(matrix)
print("\nLongest common subsequence:", aligned_seq)
