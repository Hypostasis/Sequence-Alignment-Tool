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
            print("ALIGNMENT:")
            print(self.strA)
            print(self.strB)
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
        return mi.blosum62[(query[1],query[0])]

def NeedlemanWunsch(seqA, seqB):
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
    print("Optimal value: ", H[len(seqA)][len(seqB)])
    traversal = crawler((len(seqA), len(seqB)), seqA, seqB, "","", antecedents)

    traversal.step()
    return