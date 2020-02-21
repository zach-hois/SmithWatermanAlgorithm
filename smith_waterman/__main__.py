import sys
#from .algs import sw
from .algs2 import sw
from .read_PAM import read_matrix
#the code below will allow me to take the positive and negative pairs and take the sequence
def pairs(filename):
	with open(filename) as fh:
		for line in fh:
			line = line.strip().split()
			yield line[0], line[1]
def parseFasta(filename): #This is used to take just the sequence from the fasta document
	seq = ""
	with open(filename) as fh:
		for line in fh:
			if line.startswith(">"):
				continue
			seq += line.strip()
	return seq

posMatches = [] #initializations
allFiles = []

BLOSUM50 = read_matrix("./BLOSUM50")
#print(BLOSUM50)

for file in pairs("./Pospairs.txt"): #take the file names from the pairs and get the sequence
	posMatches.append(file)
	allFiles.append(file[0])
	allFiles.append(file[1])

negMatches = []
for file in pairs("./Negpairs.txt"):
	negMatches.append(file)
	allFiles.append(file[0])
	allFiles.append(file[1])

allFiles = set(allFiles)

sequences = {file:parseFasta("./"+file) for file in allFiles}
#take out all of the sequences and return

#seq1 = "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"
#seq2 = "ANKTRELCMKSLEHAKVDTSNEARQDGIDLYKHMFENYPPLRKYFKSREEYTAEDVQNDPFFAKQGQKILLACHVLCATYDDRETFNAYTRELLDRHARDHVHMPPEVWTDFWKLFEEYLGKKTTLDEPTKQAWHEIGREFAKEINK"
posA_SCORES = []
posNormalizedScore = []
negA_SCORES = []
negNormalizedScore = []
for pos in posMatches: #this will run through all the pairs of sequences and put them in the algorithm
	#print(sequences[pos[0]])
	#print(sequences[pos[1]])
	#saving the scores so that i can call all of the positive and negative match scores
	x, y, z = sw(seq1 = sequences[pos[0]], seq2 = sequences[pos[1]],mat= BLOSUM50)
	posA_SCORES.append(z)
	posNormalizedScore.append(y)
for neg in negMatches: #this will run through all the pairs of sequences and put them in the algorithm
	#print(sequences[pos[0]])
	#print(sequences[pos[1]])
	x,y,z = sw(seq1 = sequences[neg[0]], seq2 = sequences[neg[1]],mat=BLOSUM50)
	negA_SCORES.append(z)
	negNormalizedScore.append(y)
print(posA_SCORES)
#print(posNormalizedScore)
print(negA_SCORES)
#print(negNormalizedScore)
