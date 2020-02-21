import re
#from the help of some friends, this function will
#read the matrix, see the sequence and if it is a match or mismatch, 
#and assign a score based on the matrix that is input, giving bigger penalties
#for more aggregiuous mismatches
def read_matrix(filename):

	scoring_dict = {}
	flag = False
	with open(filename) as fh:
		count = 0
		for line in fh:
			line = line.strip()
			if line.startswith("#"): continue
			elif line.startswith("A"):
				amino_acids = line.split("  ")
				flag = True
				continue

			if flag:
				values = []
				m = re.split(' ',line) #could not figure out regex for 1 and 2 spaces
				for x in m:
					try: values.append(int(x))
					except: pass

				for index,value in enumerate(values):
					scoring_dict[(amino_acids[count],amino_acids[index])] = int(value)
				count += 1

	return scoring_dict

if __name__ == "__main__":
	print(read_matrix("../BLOSUM62"))




