########### Protein Secondary Structure ############

import numpy as np
import urllib3

# RNA SEQUENCE
http = urllib3.PoolManager()
link = "https://www.ebi.ac.uk/ena/browser/api/fasta/AB063105?display=fasta"
r = http.request("GET",link)
seq = r.data.decode("utf-8") # fica em string
seq = seq.split("\n")
del seq[0]
del seq[-1]
seq = ''.join(seq)
seq = seq.replace('T','U')

#seq = "GGGGGUGUAGCUCAGUGGUAGAGCGCGUGCUUAGCAUGUACGAGGUCCCGGGUUCAAUCCCCGGCACCUCCA"
#seq = "CCCCCCCCUUAAAAAAAAGCGUUUUUUUUCCGGGGGGGG"

# PAIR DICTIONARY (possible matches)
w = {}
amin = ["G", "A", "C", "U"]

for i in range(0, len(amin)):
    for j in range(0, len(amin)):
        par = amin[i] + amin[j]
        if par == "GC" or par == "CG" or par == "AU" or par == "UA":
            w[par] = 1
        else:
            w[par] = 0


# FILL IN THE MATRIX
matrix = np.zeros((len(seq), len(seq)))
trace = np.zeros((len(seq), len(seq)))

# The main diagonal and everything below it will be zero. To implement the algorithm on element (i,j), 
# we need the score of the element below (i+1,j), the element to the left (i, j-1) and the one below 
# and to the left (i+1,j-1). To do this, we go through the columns from left to right and the rows from
# bottom to top, and only apply the algorithm if the element belongs to the upper triangle of the matrix (j>i).

for j in range(1, len(seq)):
    for i in reversed(range(0, len(seq))):

        if j > i:

            score_downleft = matrix[i+1][j-1] + w[seq[i]+seq[j]]
            score_down = matrix[i+1][j]
            score_left = matrix[i][j-1]

            score_aux = 0
            ind = 0
            for k in range(i+1, j):
                s = matrix[i][k] + matrix[k+1][j]

                if s > score_aux:
                    score_aux = s
                    ind = k
            
            values = [score_downleft, score_down, score_left, score_aux]
            matrix[i][j] = max(values)

            if score_downleft == max(values) and w[seq[i]+seq[j]] == 1:
                trace[i][j] = 1
            elif score_down == max(values):
                trace[i][j] = 2
            elif score_left == max(values):
                trace[i][j] = 3
            else:
                trace[i][j] = - ind


# TRACEBACK
list_traceback = []
def traceback(i, j):

    if j > i: # o traceback para quando se chega à diagonal

        value = int(trace[i][j])

        if value == 1 and j-i>2:
            list_traceback.append([i,j]) # adicionam-se os elementos que constituem um par na posicao i,j
            traceback(i+1, j-1)
        
        elif value == 2:
            traceback(i+1, j)
        
        elif value == 3:
            traceback (i, j-1)

        elif value < 0:
            k = - value
            traceback(i, k)
            traceback(k+1, j)

i = 0
j = len(seq) - 1      
traceback(i,j)


# STRUCTURE REPRESENTATION
structure = list('.' * len(seq))
for ind in list_traceback:
    structure[ind[0]] = '('
    structure[ind[1]] = ')'

n_basepairs = len(list_traceback)
structure = ''.join(structure)

print("\nNucleotide Sequence: ")
print(seq)

print("\nRepresentation of Secondary Structure: ")
print(structure)

print("\nNumber of Alignments Found: ")
print(n_basepairs)






        
