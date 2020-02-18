import numpy as np

match = 2 #scores taken from the suggested on wikipedia
mismatch = -1 #these decide the next step in the matrix
gap = -1 #we need to penalize a gap
seq1 = None
seq2 = None #these can be global

def sw(seq1, seq2):
	rows = len(seq1) + 1 #for the matrix, plus a lil extra for the gap
	cols = len(seq2) + 1

	#initialize the scoring matrix
	scoreMatrix, startPosition = newScoringMatrix(rows, cols)

	#call the function to find the path through the scoring matrix,
	#aligning the sequences
	alignedSeq1, alignedSeq2 = path(scoreMatrix, startPosition)

	# print the results 
	length1 = len(alignedSeq1)
 	#saving this for last lol
    return None

def newScoringMatrix(rows, cols):
	"""
	making a scoring matrix that will then be filled with the scores
	"""
	scoreMatrix = [[0 for col in range(cols)] for row in range(rows)]
	maxScore = 0
	maxPosition = None
	for i in range(1, rows):
		for j in range(1, cols):
			score = score(scoreMatrix, i, j)
			if score > maxScore:
				maxScore = score
				maxPosition = (i, j)

			score_matrix[i][j] = score #append the score to the position

	return scoreMatrix


def score(matrix, x, y):
	"""
	calculate the score for every position in the scoring matrix 
	based on the position in the table's neighbors
	"""
	similarity = match if seq1[x-1] = seq2[y-1] else mismatch #this is meta to check if the two positions in a sequence match

	diagonalScore = matrix[x-1][y-1] + similarity #we did something like this in class
	upScore = matrix[x-1][y] + gap #vertical score accounting for gaps
	leftScore = matrix[x][y-1] + gap#horizontal score accounting for gaps

    return max(0, diagonalScore, upScore, leftScore) #return the best score of those calculated, this will be the way through the matrix and match up the seqs

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
