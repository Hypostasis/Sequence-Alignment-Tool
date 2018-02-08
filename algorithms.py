from Bio.SubsMat import MatrixInfo as mi
import numpy as np
from itertools import combinations
from Bio import pairwise2
class crawler:
    def __init__(self, point, seqA, seqB, strA, strB, matrix):
        self.point = point
        self.seqA = seqA
        self.seqB = seqB
        self.strA = strA
        self.strB = strB
        self.matrix = matrix

    def step(self):

        if(self.point == (0, 0)):
            with open('output.txt', 'a') as file:
                file.write("ALIGNMENT:\n")
                file.write(self.strA+'\n')
                file.write(self.strB+'\n')
            return
        if(self.matrix[self.point[0]][self.point[1]] == 4):
            self.strA = self.seqA[self.point[0]-1] + self.strA
            self.strB = self.seqB[self.point[1] - 1] + self.strB
            self.point = (self.point[0]-1, self.point[1] - 1)
            self.step()
        if (self.matrix[self.point[0]][self.point[1]] == 2):
            #z lewej, przerwa w A
            self.strA = "-" + self.strA
            self.strB = self.seqB[self.point[1] - 1] + self.strB
            self.point = (self.point[0], self.point[1]-1)
            self.step()
        if (self.matrix[self.point[0]][self.point[1]] == 1):
            #z gory, przerwa w B
            self.strA = self.seqA[self.point[0] - 1] + self.strA
            self.strB = "-" + self.strB
            self.point = (self.point[0] -1, self.point[1])
            self.step()
        if (self.matrix[self.point[0]][self.point[1]] == 3):
            #z lewej i z gory
            newcrawler1 = crawler((self.point[0], self.point[1]-1), self.seqA, self.seqB, "-" + self.strA, self.seqB[self.point[1] - 1] + self.strB, self.matrix)
            newcrawler2 = crawler((self.point[0]-1, self.point[1]), self.seqA, self.seqB, self.seqA[self.point[0] - 1] + self.strA, "-" + self.strB, self.matrix)
            newcrawler1.step()
            newcrawler2.step()
            return
        if (self.matrix[self.point[0]][self.point[1]] == 5):
            #po skosie i z gory
            newcrawler1 = crawler((self.point[0]-1,self.point[1]-1),self.seqA, self.seqB, self.seqA[self.point[0]-1] + self.strA, self.seqB[self.point[1] - 1] + self.strB, self.matrix)
            newcrawler2 = crawler((self.point[0]-1, self.point[1]), self.seqA, self.seqB, self.seqA[self.point[0] - 1] + self.strA, "-" + self.strB, self.matrix)
            newcrawler1.step()
            newcrawler2.step()
            return
        if (self.matrix[self.point[0]][self.point[1]] == 6):
            #po skosie i z lewej
            newcrawler1 = crawler((self.point[0]-1,self.point[1]-1),self.seqA, self.seqB, self.seqA[self.point[0]-1] + self.strA, self.seqB[self.point[1] - 1] + self.strB, self.matrix)
            newcrawler2 = crawler((self.point[0], self.point[1]-1), self.seqA, self.seqB, "-" + self.strA, self.seqB[self.point[1] - 1] + self.strB, self.matrix)

            newcrawler1.step()
            newcrawler2.step()
            return
        if (self.matrix[self.point[0]][self.point[1]] == 7):
            #po skosie i z lewej i z gory
            newcrawler1 = crawler((self.point[0] - 1, self.point[1] - 1), self.seqA, self.seqB,
                                self.seqA[self.point[0] - 1] + self.strA, self.seqB[self.point[1] - 1] + self.strB,
                                self.matrix)
            newcrawler2 = crawler((self.point[0], self.point[1]-1), self.seqA, self.seqB, "-" + self.strA, self.seqB[self.point[1] - 1] + self.strB, self.matrix)
            newcrawler3 = crawler((self.point[0]-1, self.point[1]), self.seqA, self.seqB, self.seqA[self.point[0] - 1] + self.strA, "-" + self.strB, self.matrix)

            newcrawler1.step()
            newcrawler2.step()
            newcrawler3.step()
            return



def readValue(query):
    if query in mi.blosum62.keys():
        return mi.blosum62[query]
    else:
        return mi.blosum62[(query[1], query[0])]

def NeedlemanWunsch(seqA, seqB):
    with open('output.txt', 'w') as file:
        file.write("SEQUENCE ALIGNMENT TOOL OUTPUT\n")
    gap_penalty = -7
    antecedents = np.zeros((len(seqA) + 1, len(seqB) + 1), int)
    for i in range(1, len(seqB)):
        antecedents[0][i] = 2
    for i in range(1, len(seqA)):
        antecedents[i][0] = 1

    H = np.zeros((len(seqA) + 1, len(seqB) + 1), int)
    for i in range(1, len(seqB)+1):
        H[0][i] = i * gap_penalty
    for i in range(1, len(seqA)+1):
        H[i][0] = i * gap_penalty

    for i in range(1, len(seqA)+1):
        for j in range(1, len(seqB) + 1):
            optimalValue = max(H[i-1][j-1] + readValue((seqA[i-1], seqB[j-1])), H[i-1][j] + gap_penalty, H[i][j-1] + gap_penalty)
            H[i][j] = optimalValue
            #diagonal +4, right arrow +2, down arrow +1
            if optimalValue == H[i-1][j-1] + readValue((seqA[i-1], seqB[j-1])):
                #przekatna
                antecedents[i][j] = antecedents[i][j] + 4
            if optimalValue == H[i-1][j] + gap_penalty:
                #z gory
                antecedents[i][j] = antecedents[i][j] + 1
            if optimalValue == H[i][j-1] + gap_penalty:
                #z lewej
                antecedents[i][j] = antecedents[i][j] + 2
    with open('output.txt', 'a') as file:
        file.write("Optimal value: " + str(H[len(seqA)][len(seqB)]) + "\n")
    traversal = crawler((len(seqA), len(seqB)), seqA, seqB, "","", antecedents)

    traversal.step()
    return

def compute_sum_of_pairs(alignment, scoring_matrix = mi.blosum62, gap_penalty = -1):
    """Returns the sum of pairs evaluation score of the given alignment"""
    score = 0
    for i, column in enumerate(alignment.T):
        for pair in combinations(column, 2):
            if "-" in pair:
                if pair[0] == pair[1]:
                    #score of (-,-) is 0
                    continue
                score = score + gap_penalty
                continue
            score = score + readValue(pair)


    return score

def align_star(sequences, match_score = 1, mismatch_penalty = -1, gap_penalty = -1, extension_penalty = 0):

    def extend(msa_to_extend, central, aligned_seq):
        """given the sequence aligned to the center sequence, adds the aligned sequence to the msa"""
        symbols_count = len(central) - central.count("-")
        sequence_to_append = ""

        msa_pointer = 0
        central_pointer = 0
        for i in range(symbols_count + 1):
            # iterate for each string of gaps, possibly empty
            # gap counters for current central sequence of the MSA and the central sequence of the pairwise alignment
            msa_counter = 0
            central_counter = 0
            if msa_pointer < len(msa_to_extend[0]):
                while msa_to_extend[0][msa_pointer + msa_counter] == "-":
                    msa_counter += 1
                    if msa_pointer + msa_counter == len(msa_to_extend[0]):
                        break

            if central_pointer < len(central):
                while central[central_pointer + central_counter] == "-":
                    central_counter += 1
                    if central_pointer + central_counter == len(central):
                        break


            if msa_counter == central_counter:
                # same amount of gaps
                sequence_to_append += aligned_seq[central_pointer:central_pointer + central_counter + 1]
                central_pointer = central_pointer + central_counter + 1
                msa_pointer = msa_pointer + msa_counter + 1
            elif msa_counter > central_counter:
                # introduce gaps into the sequence being added to the msa
                diff = msa_counter - central_counter
                sequence_to_append += "-"*diff + aligned_seq[central_pointer:central_counter+central_pointer+1]
                msa_pointer = msa_pointer + msa_counter + 1
                central_pointer = central_pointer + central_counter + 1
            else:
                # introduce gaps into msa

                diff = central_counter - msa_counter
                for i, sequence in enumerate(msa_to_extend):
                    msa_to_extend[i] = msa_to_extend[i][:msa_pointer] + "-" * diff + msa_to_extend[i][msa_pointer:]
                sequence_to_append += aligned_seq[central_pointer:central_pointer+central_counter+1]
                central_pointer = central_pointer + central_counter + 1
                msa_pointer = msa_pointer + msa_counter + diff + 1

        msa_to_extend.append(sequence_to_append)

        return

    #find the central sequence by pairwise alignments
    matrix = np.zeros((len(sequences), len(sequences)+1))
    for i, row in enumerate(matrix):
        for j in range(len(sequences)):
            if i == j:
                continue
            matrix[i][j] = pairwise2.align.globalms(sequences[i][1], sequences[j][1],
                                                    match_score, mismatch_penalty, gap_penalty,
                                                    extension_penalty, score_only = True)
        matrix[i][len(sequences)] = np.sum(matrix[i][0:len(sequences)])
    central_sequence_score = 0
    central_sequence = -1
    for i, row in enumerate(matrix):
        if row[len(row)-1] > central_sequence_score:
            central_sequence = i
            central_sequence_score = row[len(row)-1]
    print("CENTRAL SEQUENCE", central_sequence)
    print(str(sequences[central_sequence][1]))

    #obtain pairwise alignments with the central sequence
    alignments = []
    for i, sequence in enumerate(sequences):
        if i == central_sequence:
            continue
        alignments.append(pairwise2.align.globalms(sequences[i][1],
                                                   sequences[central_sequence][1],
                                                    match_score, mismatch_penalty, gap_penalty,
                                                    extension_penalty, one_alignment_only = True)[0])
    #construct MSA iteratively
    msa = [sequences[central_sequence][1]]
    for alignment in alignments:
        central_pattern = alignment[1]
        insertion_pattern = alignment[0]
        extend(msa, central_pattern, insertion_pattern)
    #print the result to the console
    print("ALIGNMENT:")
    for x in msa:
        print(x)

    #save the result to file
    score = compute_sum_of_pairs(np.array(msa), mi.blosum62)
    print(score)
    file = open("output.txt", "w")
    for x in msa:
        file.write(str(x))
        file.write("\n")
    file.write("Sum-of-pairs score: " + str(score))
    file.close()