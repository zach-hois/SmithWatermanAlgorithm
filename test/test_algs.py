import numpy as np
import os
from smith_waterman import algs, algs2, read_PAM

BLOSUM50 = read_PAM.read_matrix("./BLOSUM50")
#print(BLOSUM50)
BLOSUM62 = read_PAM.read_matrix("./BLOSUM62")
MATIO = read_PAM.read_matrix("./MATIO")
PAM100 = read_PAM.read_matrix("./PAM100")
PAM250 = read_PAM.read_matrix("./PAM250")

def test_smithwaterman(): #test that two identical sequences return a perfect match
	#file prot-0004.fa
	seq1 = "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"
	seq2 = "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM     "

    # update this assertion
	assert algs.sw(seq1, seq2) == (142, 1.0, 284) 
	#142 out of 142 matched, 1.0 perfect match, the score returned is 284 
"""
output looks like this 
 Identities = 142/142 (100.0%), Gaps = 0/142 (0.0%)
Seq A  1     SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPE  60  
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Seq B  1     SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPE  60  

Seq A  101   ISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM  142 
             ||||||||||||||||||||||||||||||||||||||||||
Seq B  101   ISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM  142
"""

def test_roc():
	algs2.roc()
	return None






