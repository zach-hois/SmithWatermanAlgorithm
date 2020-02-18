import numpy as np

match = 2 #scores taken from the suggested on wikipedia
mismatch = -1 #these decide the next step in the matrix
gap = -1 #we need to penalize a gap
seq1 = "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"
seq2 = "ANKTRELCMKSLEHAKVDTSNEARQDGIDLYKHMFENYPPLRKYFKSREEYTAEDVQNDPFFAKQGQKILLACHVLCATYDDRETFNAYTRELLDRHARDHVHMPPEVWTDFWKLFEEYLGKKTTLDEPTKQAWHEIGREFAKEINK"
#these can be global

def sw(seq1, seq2):
	rows = len(seq1) + 1 #for the matrix, plus a lil extra for the gap
	cols = len(seq2) + 1

	#initialize the scoring matrix
	scoreMatrix, startPosition = newScoringMatrix(rows, cols)

	#call the function to find the path through the scoring matrix,
	#aligning the sequences
	alignedSeq1, alignedSeq2 = path(scoreMatrix, startPosition)

	# print the results 
	alignmentString1, idents, gaps, mismatches = alignmentString(alignedSeq1, alignedSeq2)
	length1 = len(alignedSeq1)

	#using our alignment visualization function to show what is matches/not in the sequences
	print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,length1, idents / length1, gaps, length1, gaps / length1))

	for i in range(0, length1, 60):
		seq1_slice = alignedSeq1[i:i+60]
		print('SequenceA  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
		print('             {0}'.format(alignmentString1[i:i+60]))
		seq2_slice = alignedSeq2[i:i+60]
		print('SequenceB  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
		print()
	return alignmentString1
    

def newScoringMatrix(rows, cols):
	"""
	making a scoring matrix that will then be filled with the scores
	the best alignment is the path with the highest cumulative score
	"""
	scoreMatrix = [[0 for col in range(cols)] for row in range(rows)] #initialize
	maxScore = 0
	maxPosition = None
	for i in range(1, rows): #begin looking through the matrix
		for j in range(1, cols):
			score1 = score(scoreMatrix, i, j) #calculate the score for the position
			if score1 > maxScore: #put in the score, assign to postiion
				maxScore = score1
				maxPosition = (i, j)  

			scoreMatrix[i][j] = score1 #append the score to the position

	return scoreMatrix, maxPosition


def score(matrix, x, y):
	"""
	calculate the score for every position in the scoring matrix 
	based on the position in the table's neighbors
	"""
	similarity = match if seq1[x-1] == seq2[y-1] else mismatch #this is meta to check if the two positions in a sequence match

	diagonalScore = matrix[x-1][y-1] + similarity #we did something like this in class
	upScore = matrix[x-1][y] + gap #vertical score accounting for gaps
	leftScore = matrix[x][y-1] + gap #horizontal score accounting for gaps

	maximum = max(0, diagonalScore, upScore, leftScore)

    #return the best score of those calculated, this will be the way through the matrix and match up the seqs
	return maximum 

def path(scoreMatrix, startPosition): 
	"""
	how to decide which path to take through the scoring matrix
	taking the best path through the matrix based on the score will 
	align the sequences with the best match and gap sequence/ratio

	diagonal move = match or mismatch
	up move = gap in seq1
	left move = gap in seq2
	"""

	END, DIAG, UP, LEFT = range(4) #the potential moves we can make
	alignedSeq1 = [] #initializations
	alignedSeq2 = []
	x, y = startPosition
	move = moveFunction(scoreMatrix, x, y)
	while move != END: #if there is a scheduled move according to the function
		if move == DIAG: #moving diagonally will be one up to x and y
			alignedSeq1.append(seq1[x-1])
			alignedSeq2.append(seq1[y-1])
			x-=1
			y-=1
		elif move == UP: #moving up (gap) would be just change in x
			alignedSeq1.append(seq1[x-1])
			alignedSeq2.append('-')
			x-=1
		else: #move sideways would be a change in the y here
			alignedSeq1.append('-')
			alignedSeq2.append(seq1[y-1])
			y-=1
		move = moveFunction(scoreMatrix, x, y) #keep it going in the loop

	alignedSeq1.append(seq1[x-1])
	alignedSeq2.append(seq1[y-1])

	return ''.join(reversed(alignedSeq1)), ''.join(reversed(alignedSeq2)) #return the seq

def moveFunction(scoreMatrix, x, y):
	"""
	this will be the function that allows us to move throughout the scoring matrix
	"""
	diag = scoreMatrix[x-1][y-1] #all of the possible moves based on 
	up = scoreMatrix[x-1][y] #coordinates in the scoring matrix
	left = scoreMatrix[x][y-1]

	if diag >= up and diag >= left: #in this function, a score "tie" is a diagonal move
		return 1 if diag != 0 else 0 #0 in this case is the end
	elif up > diag and up >= left: #tie here goes to the "up" move
		return 2 if up != 0 else 0
	elif left > diag and left > up:
		return 3 if left != 0 else 0 #left move or end of sequence
	else: #this is a safety measure
		raise ValueError('invalid')

def alignmentString(alignedSeq1, alignedSeq2):
	"""
	this is a visualization of how well the sequences are aligned, we'll print this
	same nucleotide is |, gaps are -, and mismatch is X
	"""
	idents, gaps, mismatches = 0,0,0 #initialization 
	alignmentString = []
	for base1, base2 in zip(alignedSeq1, alignedSeq2):
		if base1 == base2: #if the bases are the same we show the match symbol
			alignmentString.append('|')
			idents +=1
		elif '-' in (base1, base2): #if there is a gap required we put in a space
			alignmentString.append(' ')
			gaps += 1
		else:
			alignmentString.append('X') #if they do not match but no gap put an X
			mismatches += 1
	return ''.join(alignmentString), idents, gaps, mismatches


def returnMatrix(matrix):
	"""
	using this to visualize the scoring matrix for sanity
	"""
	for row in matrix:
		for col in row:
			print('{0:4}'.format(col))
			print()

def roc():
    return None
