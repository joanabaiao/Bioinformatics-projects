# NEEDLEMAN-WUNSCH ALGORITHM: Global sequence alignment

import numpy as np

seq1 = "ACGTCATCA"
seq2 = "TAGTGTCA"

matrix = np.zeros((len(seq1)+1, len(seq2)+1))
matrix_pointers = np.zeros((len(seq1)+1, len(seq2)+1))

# WEIGHTS
gap = -2  # insertion or deletion
match = 1  # match
mismatch = -1  # mismatch (mutation)

for i in range(0,len(seq1)+1):
    for j in range(0,len(seq2)+1):

        score_up = matrix[i-1][j] + gap
        score_left = matrix[i][j-1] + gap
        if seq1[i-1] == seq2[j-1]:  
            score_diagonal = matrix[i-1][j-1] + match
        else:
            score_diagonal = matrix[i-1][j-1] + mismatch

        if i != 0 or j != 0:
            if i == 0:
                matrix[i,j] = score_left
                matrix_pointers[i,j] = 2        
            elif j == 0: 
                matrix[i,j] = score_up 
                matrix_pointers[i,j] = 1
            else:
                next_value = [score_up, score_left, score_diagonal]
                matrix[i,j] = max(next_value)
                matrix_pointers[i,j] = next_value.index(max(next_value)) + 1 
                # 1 if the result comes from above, 2 from the left, and 3 from the diagonal

    
print(matrix)
print("\n", matrix_pointers)
[i,j] = np.shape(matrix)
i = i - 1
j = j - 1

align_seq1 = []
align_seq2 = []
compare_seq = []

# TRACEBACK
while i > 0 or j > 0:

    if matrix_pointers[i, j] == 1:  # Up: insertion
        align_seq1.append(" ")
        align_seq2.append(seq1[i-1])
        compare_seq.append(" ")
        i = i - 1

    elif matrix_pointers[i, j] == 2:  # Left: deletion
        align_seq1.append(seq2[j-1])
        align_seq2.append(" ")
        compare_seq.append(" ")
        j = j - 1

    elif matrix_pointers[i, j] == 3:  # Diagonal
        align_seq1.append(seq2[j-1])
        align_seq2.append(seq1[i-1])

        if seq1[i-1] == seq2[j-1]:
            compare_seq.append("|") 
        else:
            compare_seq.append("x")
        
        i = i - 1
        j = j - 1
 
align_seq1.reverse()
align_seq2.reverse()
compare_seq.reverse()
aligned_seq = " ".join(align_seq1) + "\n" + " ".join(compare_seq) + "\n" + " ".join(align_seq2)

print("\nAlignment of the sequences:")
print(aligned_seq)

