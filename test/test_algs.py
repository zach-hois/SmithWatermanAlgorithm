import numpy as np
import os
from smith_waterman import algs
"""
def test_roc():
    return None

def test_smithwaterman():
    return None
"""

def test_scoring(): #test that two identical sequences return a perfect match
	#file prot-0004.fa
	seq1 = "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"
	seq2 = "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDMXXXXX"

    # update this assertion
	assert algs.sw(seq1, seq2) == (142,142)





