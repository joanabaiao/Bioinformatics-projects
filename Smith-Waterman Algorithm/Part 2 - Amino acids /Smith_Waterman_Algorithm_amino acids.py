################ Sequence Search ###################

import urllib3
import numpy as np
import matplotlib.pyplot as plt

http = urllib3.PoolManager()

# GET SEQUENCES OF THE 2 PROTEINS
link = "http://www.uniprot.org/uniprot/?sort=&desc=&compress=no&query=cluster:(UniRef50_P51587)&fil=&force=no&preview=true&format=list"
r = http.request("GET",link)
prot_names = r.data.decode("utf-8") 
prot_names = prot_names.split("\n") 

link = "https://www.uniprot.org/uniprot/P51587.fasta"
r = http.request("GET",link)
seq1 = r.data.decode("utf-8") 
seq1 = seq1.split("\n") 
del seq1[0]
del seq1[-1]
seq1 = ''.join(seq1)

prot_nr = 3 # to choose another protein, just change the number
link = "http://www.uniprot.org/uniprot/" + prot_names[prot_nr] + ".fasta"
r = http.request("GET",link)
seq2 = r.data.decode("utf-8")
seq2 = seq2.split("\n")
del seq2[0]
del seq2[-1]
seq2 = ''.join(seq2)


# LOAD BLOSUM62 MATRIX AND FORM DICTIONARY WITH MUTATION VALUES
f = open('BLOSUM62.txt', 'r')
blosum = f.read().split("\n")
f.close()

dic_blosum = {}
aa = blosum[0].split(" ")
aa = list(filter(lambda a: a != '', aa))

for i in range(1, len(blosum)):
    values = blosum[i].split(" ")
    values = list(filter(lambda a: a != '', values))
    for j in range(1, len(values)):
        s = aa[j-1] + values[0]
        dic_blosum[s] = int(values[j])


# SMITH-WATERMAN ALGORITHM

# DEFINE SCORE AND TRACE MATRIX
score_SW = np.zeros((len(seq1)+1, len(seq2)+1))
trace_SW = np.zeros((len(seq1)+1, len(seq2)+1))

delet = -50
insert = -50

for i in range(1, len(seq1)+1):
    for j in range(1, len(seq2)+1):

        if (i == 0) and (j != 0):
            score_SW[i][j] = 0
            trace_SW[i][j] = 0

        if (j == 0) and (i != 0):
            score_SW[i][j] = 0
            trace_SW[i][j] = 0

        if (i != 0) and (j != 0):
            
            a = score_SW[i-1][j] + insert
            b = score_SW[i][j-1] + delet

            aminoacids = seq1[i-1]+seq2[j-1]
            c = score_SW[i-1][j-1] + int(dic_blosum[aminoacids])

            values = [0, a, b, c]
            score_SW[i][j] = max(values)

            if score_SW[i][j] == a:
                trace_SW[i][j] = 1

            elif score_SW[i][j] == b:
                trace_SW[i][j] = 2

            elif score_SW[i][j] == c:
                trace_SW[i][j] = 3

            elif score_SW[i][j] == 0:
                trace_SW[i][j] = 0

# print(score_SW)
# print(trace_SW)
# plt.figure()
# plt.imshow(score_SW)
# plt.colorbar(label='score')
# plt.title('Smith–Waterman Matrix')
# plt.ylabel('Sequence 2')
# plt.xlabel('Sequence 1')
# plt.show()


# TRACEBACK
indexes = []
amino_seq1 = []
amino_seq2 = []
align_start = []
align_end = []

def traceback(i, j):

    new_align = False
    align_seq1 = []
    align_seq2 = []
    compare = []
    positions = []
    verify = True

    while score_SW[i][j] != 0:
        current_position = [i,j]

        if current_position not in indexes:

            positions.append(current_position)

            if trace_SW[i][j] == 1:  # INSERÇÃO
                align_seq2.append("-")
                compare.append(" ")
                align_seq1.append(seq1[i-1])
                i = i-1

            elif trace_SW[i][j] == 2:  # DELEÇÃO
                align_seq2.append(seq2[j-1])
                compare.append(" ")
                align_seq1.append("-")
                j = j-1

            elif trace_SW[i][j] == 3:

                align_seq2.append(seq2[j-1])
                align_seq1.append(seq1[i-1])
                compare.append(":")

                i = i-1
                j = j-1

        else:
            verify = False
            break

    if verify == True and len(positions) >= 200:

        new_align = True

        for i in positions:
            indexes.append(i)

        align_seq1.reverse()
        align_seq2.reverse()
        compare.reverse()

        align_seq1 = ''.join(align_seq1)
        align_seq2 = ''.join(align_seq2)

        amino_seq1.append(align_seq1)
        amino_seq2.append(align_seq2)

    return new_align, positions


new_align = False
i = len(seq1)
j = len(seq2)
cmax = len(seq2)

while (i > 0) or (j > 0):

    if (score_SW[i][j] > 0):

        if (i == len(seq1) and j == len(seq2)):
            if (score_SW[i][j] > max(score_SW[i-1][j-1], score_SW[i-1][j], score_SW[i][j-1])):
                new_align, positions = traceback(i, j)

        elif (i == len(seq1) and j != len(seq2)):
            if (score_SW[i][j] > max(score_SW[i-1][j-1], score_SW[i-1][j], score_SW[i][j-1], score_SW[i][j+1])):
                new_align, positions = traceback(i, j)

        elif (i != len(seq1) and j == len(seq2)):
            if (score_SW[i][j] > max(score_SW[i-1][j-1], score_SW[i-1][j], score_SW[i][j-1], score_SW[i+1][j])):
                new_align, positions = traceback(i, j)

        elif (i != len(seq1) and j != len(seq2)):
            if (score_SW[i][j] > max(score_SW[i-1][j-1], score_SW[i-1][j], score_SW[i][j-1], score_SW[i+1][j], score_SW[i+1][j+1], score_SW[i][j+1])):
                new_align, positions = traceback(i, j)

    if new_align == True:

        [i, j] = positions[-1]
        [rmax, cmax] = positions[-1]

        align_start.append(positions[0])
        align_end.append(positions[-1])
        
        new_align = False

    else:
        if (j == 0):
            i = i-1
            j = cmax

        else:
            j = j-1


# IMAGE OF THE ALIGNMENTS
matrix_align = np.zeros((len(seq1)+1, len(seq2)+1), dtype=np.bool_)

for coord in indexes:
    c1 = coord[0]
    c2 = coord[1]
    matrix_align[c1, c2] = True

plt.matshow(matrix_align, fignum=1, cmap='Greys')
plt.imsave('fig_part2.png', matrix_align, cmap='Greys')


print("\n------------------------------------------")
print("               ALIGNMENTS                ")
print("------------------------------------------\n")

for i in range(0, len(align_start)):
    print('ALIGNMENT', i+1, "\n")
    print("Sequence of protein P51587:")
    print(amino_seq1[i])
    print("\nProtein sequence", prot_names[prot_nr],":")
    print(amino_seq2[i], "\n\n")