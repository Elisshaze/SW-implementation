"""
Project for the Algorithms for Bioinformatics course
Author: Elisa Pettinà
Contact: elisa.pettina@studenti.unitn.it
Repository: https://github.com/Elisshaze/SW-implementation
"""
#Importing needed libriaries
import sys, getopt
import numpy as np
import itertools

def findall(element, matrix):
	result = []
	for i in range(len(matrix)):
		for j in range(len(matrix[i])):
			if matrix[i][j] == element:
				result.append((i, j))
	return result

def return_trace_score(list_moves, choice):
	#For each move, a code is assigned
	#match, insert, delete, 0
	#	0	 1	  2
	supp=[]
	trace=0
	for i, x in enumerate(list_moves):
		if choice == x:
			supp.append(i)

	if len(supp) == 3 :
		trace = 6 #all three
	elif len(supp) == 1 :
		trace = supp[0] #just one
	elif len(supp) == 2:
		if (supp[0] + supp[1]) == 1 :
			trace = 3 #match and insert
		elif (supp[0] + supp[1]) == 2 :
			trace = 4 #match and deletion
		elif (supp[0] + supp[1]) == 3 :
			trace = 5 #delete and insert

	return trace

def calculate_matrix(seq1, seq2):
	#Compute the scoring matrix
	M = np.zeros((len(seq1) + 1, len(seq2) + 1), dtype=int)
	#Compute the tracing matrix, will be used for the backtrack
	M_trace = np.zeros((len(seq1) + 1, len(seq2) + 1), dtype=int)

	for i, j in itertools.product(range(1, M.shape[0]), range(1, M.shape[1])):
		match = M[i - 1, j - 1] + (score if seq1[i - 1] == seq2[j - 1] else - score)
		insert = M[i, j - 1] - gap
		delete = M[i - 1, j] - gap
		choice = max(match, insert,  delete, 0)
		M[i, j] = choice
		#M_trace is used for backtrack
		M_trace [i, j] = return_trace_score([match, insert, delete, 0], choice)

	return M, M_trace

def findall(element, matrix):
	#Find indexes of each recurrence of an element int he matrix
	result = []
	for i in range(len(matrix)):
		for j in range(len(matrix[i])):
			if matrix[i][j] == element:
				result.append((i, j))
	return result

def find_other_max(M, M_rev, next_max, many=3):
	#Function to find all the other maxes
	M_flat = M.flatten()
	M_flat.sort()

	total = []

	for i in range(0, many) :
		second_max = M_flat[(-next_max)-1] if len(M_flat) > 1 else M_flat[0]
		indexes = findall(second_max, M_rev)
		next_max = (next_max)+len(indexes)

		total.append((second_max, indexes))

	print(total)

	return(total)

def wrapper_restore(M, M_trace, seq2):
	#Wrapper to pass the starting points to the recursive algorithm
	b_="" #List that will be popoulated in the backtrack

	M_rev = np.flipud(np.fliplr(M))
	M_trace_rev = np.flipud(np.fliplr(M_trace))

	i, j = np.unravel_index(M_rev.argmax(), M_rev.shape)

	start_points = findall(M_rev[i][j], M_rev)
	if verbose == 1 :
		print ("Starting points:", start_points) # Will return list of tuples

	for x in start_points :
		print("Score: ", M_rev[i][j], "Match score:", score, "Gap penalty:", gap)

		if M_rev[i][j] == 0 :
			print("\033[93mNo possible alignment for starting point with value", 0, '\033[0m')
			break

		restore_sequence(M_rev, M_trace_rev, seq2, b_, x[0], x[1])
		#seq.append(M_rev[i][j])
		#print('\n')

	#In case of wanting the n best score alignments
	if other_max > 0 :
		other_max_sequences = find_other_max(M, M_rev, len(start_points), other_max)
		for x in other_max_sequences :
			start = x[0]
			start_points=x[1]
			for x in start_points :
				print("Score: ", start, "Match score:", score, "Gap penalty:", gap)
				restore_sequence(M_rev, M_trace_rev, seq2, b_, x[0], x[1])
				#print('\n')

def make_move(index, b_, b, current_max, old_max,j, m_c, mm_c, d_c, i_c):
	#Update string and parameters based on the 
	if index == 0 :
		if current_max == old_max + 3 :
			b_ = "X" + b_
			mm_c = mm_c + 1
		else:
			b_= b[len(b)-j-1] + b_
			m_c = m_c + 1

	elif index  == 1 :
		b_ = "+" + b_
		i_c = i_c + 1

	elif index == 2 :
		b_ = "-" + b_
		d_c = d_c + 1

	return b_, m_c, mm_c, d_c, i_c

def restore_sequence(M_rev, M_trace_rev, b, b_, i, j, m_c=0, mm_c=0, i_c=0, d_c=0):
	l=[]
	#Base case
	if M_rev[i, j] == 0:
		return[b_, (len(b)-j), m_c, mm_c, i_c, d_c]

	old_max = M_rev[i][j]

	#Now, the possible choices are based on the origin of the best result
	# 0 -> from match
	# 1 -> from insertion
	# 2 -> from deletion
	# 3 -> from match and deletion
	# 4 -> from match and insertion
	# 5 -> from insertion and deletion
	# 6 -> from all

	if M_trace_rev[i][j] == 0:
		#match, diagonal
		current_max = M_rev[i+1][j+1]
		b_, m_c, mm_c, d_c, i_c = make_move(0,b_,b,current_max,old_max,j,m_c,mm_c,d_c,i_c)
		if verbose == 1 :
			print("a match, old max and current max:", old_max, current_max)
			print(b_)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i+1, j+1, m_c, mm_c, i_c, d_c)
		l.append(seq)

	if M_trace_rev[i][j] == 1:
		#insert
		current_max = M_rev[i][j+1]
		b_, m_c, mm_c, d_c, i_c = make_move(1,b_,b,current_max,old_max,j,m_c,mm_c,d_c,i_c)
		
		if verbose == 1 :
			print("a insertion, old max and current max:", old_max, current_max)
			print(b_)

		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i, j+1, m_c, mm_c, i_c, d_c)
		l.append(seq)


	if M_trace_rev[i][j] == 2:
		#deletion
		current_max = M_rev[i+1][j]
		b_, m_c, mm_c, d_c, i_c = make_move(2,b_,b,current_max,old_max,j,m_c,mm_c,d_c,i_c)
		
		if verbose == 1 :
			print("a deletion, old max and current max:", old_max, current_max)
			print(b_)
		
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i+1, j, m_c, mm_c, i_c, d_c)
		l.append(seq)

	#FROM MATCH AND DELETE
	if M_trace_rev[i][j] == 4:
		supp = b_
		supp_m_c, supp_mm_c, supp_d_c, supp_i_c = m_c, mm_c, d_c, i_c
		if verbose == 1 :
			print("in 4 (match and delete), old max:", old_max)

		#going up\down, deletion
		current_max = M_rev[i+1][j]
		b_, m_c, mm_c, d_c, i_c = make_move(2,b_,b,M_rev[i+1][j],old_max,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_,i+1,j, m_c, mm_c, i_c, d_c)
		l.append(seq)

		#going diagonal
		current_max = M_rev[i+1][j+1]
		b_ = supp
		m_c, mm_c, d_c, i_c = supp_m_c, supp_mm_c, supp_d_c, supp_i_c
		b_, m_c, mm_c, d_c, i_c = make_move(0,b_,b,M_rev[i+1][j+1],old_max,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i+1,j+1, m_c, mm_c, i_c, d_c)
		l.append(seq)

	#FROM MATCH AND INSERTION
	if M_trace_rev[i][j] == 3:
		supp = b_
		supp_m_c, supp_mm_c, supp_d_c, supp_i_c = m_c, mm_c, d_c, i_c

		if verbose == 1 :
			print("in 3 (match and insertion), old max:", old_max)

		#going left
		current_max = M_rev[i][j+1]
		b_, m_c, mm_c, d_c, i_c = make_move(1,b_,b,M_rev[i][j+1],old_max,j-1,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i, j+1, m_c, mm_c, i_c, d_c)
		l.append(seq)

		#going diagonal
		current_max = M_rev[i+1][j+1]
		b_ = supp
		m_c, mm_c, d_c, i_c = supp_m_c, supp_mm_c, supp_d_c, supp_i_c
		b_, m_c, mm_c, d_c, i_c = make_move(0,b_,b,M_rev[i+1][j+1],old_max,j+1,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i+1,j+1, m_c, mm_c, i_c, d_c)
		l.append(seq)

	#FROM DELETION AND INSERTION
	if M_trace_rev[i][j] == 5:
		supp = b_
		supp_m_c, supp_mm_c, supp_d_c, supp_i_c = m_c, mm_c, d_c, i_c

		if verbose == 1 :
			print("in 5 (insertion and deletion), old max:", old_max)

		#going left
		current_max = M_rev[i][j+1]
		b_, m_c, mm_c, d_c, i_c = make_move(1,b_,b,M_rev[i][j+1],old_max,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i, j+1, m_c, mm_c, i_c, d_c)
		l.append(seq)

		#going up/down
		current_max = M_rev[i+1][j]
		b_ = supp
		m_c, mm_c, d_c, i_c = supp_m_c, supp_mm_c, supp_d_c, supp_i_c
		b_, m_c, mm_c, d_c, i_c = make_move(2,b_,b,M_rev[i+1][j],old_max,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i+1, j, m_c, mm_c, i_c, d_c)
		l.append(seq)
	#FROM MATCH, INSERTION, DELETION
	if M_trace_rev[i][j] == 6:
		if verbose == 1 :
			print("in 6 (match,deletion and insertion), old max:", old_max)
		supp = b_
		supp_m_c, supp_mm_c, supp_d_c, supp_i_c = m_c, mm_c, d_c, i_c

		#going left
		current_max = M_rev[i][j+1]
		b_, m_c, mm_c, d_c, i_c = make_move(1,b_,b,M_rev[i][j+1],old_max,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i, j+1, m_c, mm_c, i_c, d_c)
		l.append(seq)

		#going up/down
		current_max = M_rev[i+1][j]
		b_ = supp
		m_c, mm_c, d_c, i_c = supp_m_c, supp_mm_c, supp_d_c, supp_i_c
		b_, m_c, mm_c, d_c, i_c = make_move(2,b_,b,M_rev[i+1][j],old_max,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i+1, j, m_c, mm_c, i_c, d_c)
		l.append(seq)

		#going diagonal
		current_max = M_rev[i+1][j+1]
		b_ = supp
		m_c, mm_c, d_c, i_c = supp_m_c, supp_mm_c, supp_d_c, supp_i_c
		b_, m_c, mm_c, d_c, i_c = make_move(0,b_,b,M_rev[i+1][j+1],old_max,j+1,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i+1,j+1, m_c, mm_c, i_c, d_c)
		l.append(seq)

	#Printing the results in stdout
	if current_max == 0:
		for x in l:
			results.append(l)
			print("SEQUENCE: ", x[0])
			if long_output == 1 :
				print("starting_position: ", x[1])
				print("matches: ", x[2])
				print("mismatches: ", x[3])
				print("insertions: ", x[4])
				print("deletions: ", x[5])
			print("------------------------------------------------------")
			print('\n')



def smith_watermann(seq1, seq2) :
	#The smith waterman is composed by a dynamic programming
	#part and a backtrack part
	M, M_trace = calculate_matrix(seq1, seq2)
	if verbose == 1 :
		print('Scoring matrix:\n', M, 'Tracing matrix\n', M_trace)

	print("\n\033[94mRESULTS:\033[0m \n")
	wrapper_restore(M, M_trace, seq2)

def args_options():
	score=3
	gap=2
	other_max=0
	long_output = 0
	verbose = 0
	argumentList = sys.argv[1:]
	case=0

	#Options
	options = "hlvm:g:b:c"
	#Long options
	long_options = ["help", "long_output", "verbose", "match", "gap", "best", "case"]

	try:
		arguments, values = getopt.getopt(argumentList, options, long_options)
		# checking each argument
		for currentArgument, currentValue in arguments:
			if currentArgument in ("-l", "--long_ouput"):
				long_output = 1

			if currentArgument in ("-v", "--verbose"):
				verbose = 1

			if currentArgument in ("-m", "--match"):
				score = int(currentValue)
				print("Match score / mismatch penalty:", score)

			if currentArgument in ("-g", "--gap"):
				gap = int(currentValue)
				print("Gap penalty:", gap)

			if currentArgument in ("-b", "--best"):
				other_max = int(currentValue)

			if currentArgument in ("-c", "--case"):
				case=1
			if currentArgument in ("-h", "--help"):
				with open('help.txt', 'r') as f:
					print(f.read())
				exit()

	except getopt.error as err:
		#Output error, and return with an error code
		print (str(err))

	return score, gap, other_max, long_output, verbose, case


def user_input():
	seq1 = input('Type reference sequence:\n')
	if seq1:
		print("\033[92mOK!\033[0m")
	else :
		print("insert a valid sequence")
		exit()
	seq2 = input('Type matching sequence:\n')
	if seq2:
		print("\033[92mOK!\033[0m")
	else :
		print("insert a valid sequence")
		exit()

	return(seq1, seq2)


#The matching and gap scores and how many best alignment are defined by options
#Options' meanings are accessible thour -h or in the help.txt file
score, gap, other_max, long_output, verbose, case = args_options()
#The user inputs the sequences when asked
seq1, seq2 = user_input()

#Array for the results
results = []

if case==1 :
	seq1 = seq1.upper()
	seq2 = seq2.upper()
	print(seq1, seq2)

smith_watermann(seq1, seq2)

#if you want to print the results outside the recursion:
#print("all results", results)
