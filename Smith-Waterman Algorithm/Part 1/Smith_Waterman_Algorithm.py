################ Sequence Search ###################

import numpy as np
import matplotlib.pyplot as plt

# LOAD DATA SEQUENCES
dna = open("DNA_sequence.txt", "r")
dna = dna.read()
dna = dna.split("\n")
del dna[0]
del dna[len(dna)-2:len(dna)]
dna = ''.join(dna)

mrna = open("mRNA_sequence.txt", "r")
mrna = mrna.read()
mrna = mrna.split("\n")
del mrna[0]
del mrna[len(mrna)-2:len(mrna)]
mrna = ''.join(mrna)


# SMITH-WATERMAN ALGORITHM

# DEFINE SCORE AND TRACE MATRIX
score_SW = np.zeros((len(dna)+1, len(mrna)+1))
trace_SW = np.zeros((len(dna)+1, len(mrna)+1))

delet = -50
insert = -50
match = +5
mismatch = -5

for i in range(1, len(dna)+1):
    for j in range(1, len(mrna)+1):

        if (i == 0) and (j != 0):
            score_SW[i][j] = 0
            trace_SW[i][j] = 0

        if (j == 0) and (i != 0):
            score_SW[i][j] = 0
            trace_SW[i][j] = 0

        if (i != 0) and (j != 0):
            a = score_SW[i-1][j] + insert
            b = score_SW[i][j-1] + delet

            if dna[i-1] == mrna[j-1]:
                c = score_SW[i-1][j-1] + match
            else:
                c = score_SW[i-1][j-1] + mismatch

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

print(score_SW)
# print(trace_SW)
# plt.figure()
# plt.imshow(score_SW)
# plt.colorbar(label='score')
# plt.title('Smith–Waterman Matrix')
# plt.ylabel('mRNA Sequence')
# plt.xlabel('DNA Sequence')
# plt.show()


# TRACEBACK
indexes = []
exons_mrna = []
exons_dna = []
exons_start = []
exons_end = []

def traceback(i, j):

    new_exon = False
    exon_seq = []
    dna_seq = []
    compare = []
    positions = []
    verify = True

    while score_SW[i][j] != 0:
        current_position = [i,j]

        if current_position not in indexes:

            positions.append(current_position)

            if trace_SW[i][j] == 1:  # INSERÇÃO
                exon_seq.append("-")
                compare.append(" ")
                dna_seq.append(dna[i-1])
                i = i-1

            elif trace_SW[i][j] == 2:  # DELEÇÃO
                exon_seq.append(mrna[j-1])
                compare.append(" ")
                dna_seq.append("-")
                j = j-1

            elif trace_SW[i][j] == 3:

                if dna[i-1] == mrna[j-1]:
                    exon_seq.append(mrna[j-1])
                    dna_seq.append(dna[i-1])
                    compare.append("|")

                else:
                    exon_seq.append(mrna[j-1])
                    dna_seq.append(dna[i-1])
                    compare.append("x")

                i = i-1
                j = j-1

        else:
            verify = False
            break

    if verify == True and len(positions) >= 200:

        new_exon = True

        for i in positions:
            indexes.append(i)

        exon_seq.reverse()
        compare.reverse()
        dna_seq.reverse()

        exon_seq = ''.join(exon_seq)
        dna_seq = ''.join(dna_seq)

        exons_mrna.append(exon_seq)
        exons_dna.append(dna_seq)

    return new_exon, positions


new_exon = False
i = len(dna)
j = len(mrna)

while (i > 0) or (j > 0):

    if (score_SW[i][j] > 0):

        if (i == len(dna) and j == len(mrna)):
            if (score_SW[i][j] > max(score_SW[i-1][j-1], score_SW[i-1][j], score_SW[i][j-1])):
                new_exon, positions = traceback(i, j)

        elif (i == len(dna) and j != len(mrna)):
            if (score_SW[i][j] > max(score_SW[i-1][j-1], score_SW[i-1][j], score_SW[i][j-1], score_SW[i][j+1])):
                new_exon, positions = traceback(i, j)

        elif (i != len(dna) and j == len(mrna)):
            if (score_SW[i][j] > max(score_SW[i-1][j-1], score_SW[i-1][j], score_SW[i][j-1], score_SW[i+1][j])):
                new_exon, positions = traceback(i, j)

        elif (i != len(dna) and j != len(mrna)):
            if (score_SW[i][j] > max(score_SW[i-1][j-1], score_SW[i-1][j], score_SW[i][j-1], score_SW[i+1][j], score_SW[i+1][j+1], score_SW[i][j+1])):
                new_exon, positions = traceback(i, j)

    if new_exon == True:

        [i, j] = positions[-1]
        [rmax, cmax] = positions[-1]

        exons_start.append(positions[0])
        exons_end.append(positions[-1])
        
        new_exon = False

    else:
        if (j == 0):
            i = i-1
            j = cmax

        else:
            j = j-1


# IMAGE OF THE ALIGNMENTS
matrix_align = np.zeros((len(dna)+1, len(mrna)+1), dtype=np.bool_)

for coord in indexes:
    c1 = coord[0]
    c2 = coord[1]
    matrix_align[c1, c2] = True

plt.matshow(matrix_align, fignum=1, cmap='Greys')
plt.imsave('fig_part1.png', matrix_align, cmap='Greys')


# INTRONS AND EXONS
exons_mrna.reverse()
exons_dna.reverse()
exons_start.reverse()
exons_end.reverse()

introns = []
introns_start = []
introns_end = []
for i in range(0,len(exons_start)):
   
    if (i == 0):
        start = 0
        end = exons_end[0][0] - 1
    
    else:
        start = exons_start[i-1][0] + 1
        end = exons_end[i][0] - 1
        
    intron = dna[start:end]
    introns.append(intron)
    introns_start.append(start)
    introns_end.append(end)

print("\n------------------------------------------")
print("          LIST OF EXONS IN RNA            ")
print("------------------------------------------\n")

for i in range(0, len(exons_start)):
    print('EXON', i+1, "\n")
    print("Start Position:", exons_end[i][1])
    print("End Position:", exons_start[i][1])
    print("Sequence: ")
    print(exons_mrna[i], "\n\n")


print("\n------------------------------------------")
print("          LIST OF EXONS IN DNA            ")
print("------------------------------------------\n")

for i in range(0, len(exons_start)):
    print('EXON', i+1, "\n")
    print("Start Position:", exons_end[i][0])
    print("End Position:", exons_start[i][0])
    print("Sequence: ")
    print(exons_dna[i], "\n\n")


print("\n------------------------------------------")
print("         LIST OF INTRONS IN DNA           ")
print("------------------------------------------\n")

for i in range(0, len(introns)):
    print('INTRON', i+1, "\n")
    print("Start Position:", introns_start[i])
    print("End Position:", introns_end[i])
    print("Sequence: ")
    print(introns[i], "\n\n")
