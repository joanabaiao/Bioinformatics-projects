# EXAME MODELO

# Exercicio 3: A matriz M contém o alinhamento de duas sequências.
# Proponha um código que possibilite determinar o alinhamento de maior dimensão.

import numpy as np

seq1 = "ACGTCATCA"
seq2 = "TAGTGTCA"

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
            matrix[i][j] = max(matrix[i-1,j], matrix[i,j-1],  matrix[i-1,j-1])   
                     
print(matrix)
print()

def traceback(matrix):

    i = matrix.shape[0]-1
    j = matrix.shape[1]-1
    align_seq1 = []
    align_seq2 = []
    align = [] 

    c = 0
    max_align = 0

    while i != 0 and j != 0:

        score_up = matrix[i-1,j]
        score_left = matrix[i, j-1]
        score_diagonal = matrix[i-1,j-1]

        if seq1[i-1] == seq2[j-1]:
            align_seq1.append(seq1[i-1])
            align_seq2.append(seq2[j-1])
            align.append("|")
            c = c + 1
            
            i = i-1
            j = j-1

        else:
            c = 0
            value = max(score_up, score_left, score_diagonal)

            if value == score_diagonal:
                align_seq1.append(seq2[j-1])
                align_seq2.append(seq1[i-1])
                align.append("x")
                i = i-1
                j = j-1          

            elif value == score_up:
                align_seq1.append("-")
                align_seq2.append(seq1[i-1])
                align.append(" ")                
                i = i-1

            elif value == score_left:
                align_seq1.append(seq2[j-1])
                align_seq2.append("-")
                align.append(" ")
                j = j-1
        
        if c > max_align:
            max_align = c

    align_seq1.reverse()
    align_seq2.reverse()
    align.reverse()
    
    aligned_seq = " ".join(align_seq1) + "\n" + " ".join(align) + "\n" + " ".join(align_seq2)
    
    return max_align, aligned_seq


max_align, aligned_seq = traceback(matrix)
print(aligned_seq)
print("\nDimension of the longest correct alignment:", max_align)
