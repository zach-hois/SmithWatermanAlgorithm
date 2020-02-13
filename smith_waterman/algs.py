import numpy as np

match = 2 #scores taken from the suggested on wikipedia
mismatch = -1 #these decide the next step in the matrix
gap = -1
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
	


    return None

def roc():
    return None
