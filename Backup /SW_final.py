from operator import index
import sys, getopt
import numpy as np
import itertools


def user_input():
	#Function to take in input the sequences
	seq1 = input('Type reference sequence:\n')
	seq2 = input('Type matching sequence:\n')
	
	return(seq1, seq2)

def args_options():
	#functions defininf the possible options
	score=3
	gap=2
	argumentList = sys.argv[1:]
	other_max=0
	print(sys.argv)
 
	# Options
	options = "hlvm:g:b:"
	# Long options
	long_options = ["help", "long_output", "verbose", "match", "gap", "best"]
 
	try:
		# Parsing argument
		arguments, values = getopt.getopt(argumentList, options, long_options)
		# checking each argument
		for currentArgument, currentValue in arguments:

			if currentArgument in ("-l", "--long_ouput"):
				print ("Displaying:", sys.argv[0])
			
			if currentArgument in ("-v", "--verbose"):
				print ("Displaying:", sys.argv[0])
			
			if currentArgument in ("-m", "--match"):
				score = currentValue
		
			if currentArgument in ("-g", "--gap"):
				gap = currentValue

			if currentArgument in ("-b", "--best"):
				other_max = currentValue

	except getopt.error as err:
		# output error, and return with an error code
		print (str(err))

	return (score, gap, other_max)

def smith_watermann(seq1, seq2) :
	#The smith waterman is composed by a dynamic programming
	#part and a backtrack part
	M, M_trace = calculate_matrix(seq1, seq2)
	print('\n', M, '\n', M_trace)
	
	wrapper_restore(M, M_trace, seq2)

##################################################################
# Dynamic programming block #
##################################################################
def return_trace_score(list_moves, choice):
	#A 'trace' is assigned to each cell, based on its origin
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
			trace = 3 #match and deletion
		elif (supp[0] + supp[1]) == 2 :
			trace = 4 #match and insert
		elif (supp[0] + supp[1]) == 3 :
			trace = 5 #delete and insert
	return trace

def calculate_matrix(seq1, seq2):
	#A matrix is generated with the actual scores, the other one with 
	# the traces, which will be used during the backtrack part
	M = np.zeros((len(seq1) + 1, len(seq2) + 1), dtype=np.int)
	M_trace = np.zeros((len(seq1) + 1, len(seq2) + 1), dtype=np.int)  

	for i, j in itertools.product(range(1, M.shape[0]), range(1, M.shape[1])):
		match = M[i - 1, j - 1] + (score if seq1[i - 1] == seq2[j - 1] else - score)
		delete = M[i - 1, j] - gap
		insert = M[i, j - 1] - gap
		choice = max(match, delete, insert, 0)
		M[i, j] = choice
		#M_trace is used for backtrack
		M_trace [i, j] = return_trace_score([match, delete, insert, 0], choice)

	return M, M_trace
	
##################################################################
# Backtracking block #
##################################################################
def wrapper_restore(M, M_trace, seq2):
	#This is a wrapper for the recursive function, selects the starting points
	b_="" #list to populate in the backtrack
	
	#The flipping is perfomed to exploit numpy's methods
	M_rev = np.flipud(np.fliplr(M))
	M_trace_rev = np.flipud(np.fliplr(M_trace))

	i, j = np.unravel_index(M_rev.argmax(), M_rev.shape) 
	print ("max indexes:", i,j, "current max:", M_rev[i][j] )

	#Start points finds all the maxes
	start_points = findall(M_rev[i][j], M_rev)

	print ("how many max:", len(start_points))
	print (start_points) # Will return list of tuples
	
	for x in start_points :
		restore_sequence(M_rev, M_trace_rev, seq2, b_, x[0], x[1])
		print("\n \n \n \n end")
		

	'''
	if other_max > 0 :
		other_max_sequences = find_other_max(M, M_rev, len(start_points), other_max)
		for x in other_max_sequences :
			current_max = x[0]
			start_points=x[1]
			for x in start_points :
				
				seq = restore_sequence(M_rev, seq2, b_, x[0], x[1])
				total.append(current_max)
				total.append(seq)
	'''
def find_other_max(M, M_rev, next_max, many=3):
	total = []
	#Looks for the best x scores and all of their occurences
	M_flat = M.flatten()
	M_flat.sort()

	for i in range(0, many) :
		second_max = M_flat[(-next_max)-1] if len(M_flat) > 1 else M_flat[0]
		indexes = findall(second_max, M_rev)
		next_max = (next_max)+len(indexes)
		total.append((second_max, indexes))
	
	print(total)
	return(total)	

def findall(element, matrix):
	#simply returns all of the occurences in a matrix
    result = []
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] == element:
                result.append((i, j))
    return result

def restore_sequence(M_rev, M_trace_rev, b, b_, old_i, old_j, m_c=0, mm_c=0, i_c=0, d_c=0):
	print (b_)
	l=[]
	if M_rev[old_i, old_j] == 0:
		return[(b_, (len(b)-old_j), m_c, mm_c, i_c, d_c)]

	#Possible moves during backtrack
	match = M_rev[old_i+1][old_j+1]
	insert= M_rev[old_i][old_j+1]
	delete= M_rev[old_i+1][old_j]

	moves = [match, delete, insert]
	old_max = M_rev[old_i][old_j] #Previous max from precend call
	current_max = max(moves) #Current max 
	index_move = np.argmax(moves) #Which move was chosen

	if index_move == 0 : #a match
		i = old_i +1
		j = old_j +1
		print("a match")
	elif index_move == 1 : #a insert
		i = old_i
		j = old_j +1
		print("a insert")
	elif index_move == 2 : #a deletion
		i = old_i +1
		j = old_j
		print("a deletion")

	#Now, the possible choices are based on the origin of the best result
	# 0 -> from match
	# 1 -> from insertion
	# 2 -> from deletion
	# 3 -> from match and deletion
	# 4 -> from match and insertion
	# 5 -> from insertion and deletion
	# 6 -> from all

	if M_trace_rev[old_i][old_j] == 0 :
		print ("match for ", current_max)
		b_, m_c, mm_c, d_c, i_c = make_move(index_move,b_,b,current_max,old_max,i,j,m_c,mm_c,d_c,i_c) 
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i, j, m_c, mm_c, i_c, d_c)
		l.append(seq)

	elif M_trace_rev[old_i][old_j] == 1 :
		print ("insert for ", current_max)
		b_, m_c, mm_c, d_c, i_c = make_move(index_move,b_,b,current_max,old_max,i,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i, j, m_c, mm_c, i_c, d_c)
		l.append(seq)

	elif M_trace_rev[old_i][old_j] == 2 :
		print ("delete for ", current_max)
		b_, m_c, mm_c, d_c, i_c = make_move(index_move,b_,b,current_max,old_max,i,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i, j, m_c, mm_c, i_c, d_c)
		l.append(seq)

	#Cases with two options (3, 4, 5)
	if M_trace_rev[i][j] == 3 :
		#going up
		supp = b_
		print("going up")
		b_, m_c, mm_c, d_c, i_c = make_move(2,b_,b,current_max,old_max,i,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i+1, j, m_c, mm_c, i_c, d_c)
		l.append(seq)

		#going diagonal
		b_=supp
		print("going diagonal")
		b_, m_c, mm_c, d_c, i_c = make_move(0,b_,b,current_max,old_max,i,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i+1, j+1, m_c, mm_c, i_c, d_c)
		l.append(seq)
		
	elif M_trace_rev[i][j] == 4 :
		#going right
		supp = b_
		print("going right")
		b_, m_c, mm_c, d_c, i_c = make_move(1,b_,b,current_max,old_max,i,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i, j+1, m_c, mm_c, i_c, d_c)
		l.append(seq)
		
		#going diagonal
		b_=supp
		print("going diagonal")
		b_, m_c, mm_c, d_c, i_c = make_move(0,b_,b,current_max,old_max,i,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i+1, j+1, m_c, mm_c, i_c, d_c)
		l.append(seq)

	elif M_trace_rev[i][j] == 5 :
		#going right
		supp = b_
		print("going right")
		b_, m_c, mm_c, d_c, i_c = make_move(1,b_,b,current_max,old_max,i,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i, j+1, m_c, mm_c, i_c, d_c)
		l.append(seq)

		#going up
		supp = b_
		print("going up")
		b_, m_c, mm_c, d_c, i_c = make_move(2,b_,b,current_max,old_max,i,j,m_c,mm_c,d_c,i_c)
		seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i+1, j, m_c, mm_c, i_c, d_c)
		l.append(seq)

	seq=restore_sequence(M_rev,  M_trace_rev, b, b_, i, j, m_c, mm_c, i_c, d_c)
	l.append(seq)
	print (b_)
	if l[0] is not None :
		print("LISTAAAAAAAAAAAAAAA", l)

	
def make_move(index, b_, b, current_max, old_max, i, j, m_c, mm_c, d_c, i_c):
	if index == 0 :
		print("in make match")
		if current_max == old_max + 3 :
			b_ = "X" + b_
			mm_c = mm_c + 1
		else:
			b_= b[len(b)-j] + b_
			m_c = m_c + 1 

	elif index  == 1 :
		print("in make insertion")
		b_ = "+" + b_
		i_c = i_c + 1

	elif index == 2 :
		print("in make deletion")
		b_ = "-" + b_
		d_c = d_c + 1

	return b_, m_c, mm_c, d_c, i_c 

def update_counter():
	print("counters")



#MAIN
score, gap, other_max = args_options()
#seq1, seq2 = user_input()
seq1="GGTTGACTACCCGTACGTTT"
seq2="ACGGTGATGCCCGTTT"
print (len(seq1), len(seq2))
print (score, gap)

smith_watermann(seq1, seq2)



'''

match = M_rev[old_i+1][old_j+1]
	insert= M_rev[old_i][old_j+1]
	delete= M_rev[old_i+1][old_j]
	moves = [match, insert, delete]

	old_max = M_rev[old_i][old_j]

	

	current_max = max(moves)
	index_move = np.argmax(moves)

	print("CONTROLLO: MAX ", old_max, current_max)
	
	if index_move == 0 : #a match
		i = old_i +1
		j = old_j +1
	elif index_move == 1 : #a insert
		i = old_i
		j = old_j +1
	elif index_move == 2 : #a deletion
		i = old_i +1
		j = old_j

		'''