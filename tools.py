#! /usr/bin/python3 

import math
import os
import subprocess as sp


if __name__ == '__main__':
	print("This library is not meant to be run on its own!")
	exit()


#########################################
# Functions to facilitate spin flipping #
#########################################

def binomial(x, y):
    try:
        binom = math.factorial(x) // math.factorial(y) // math.factorial(x - y)
    except ValueError:
        binom = 0
    return binom


def getCombinations(digs,places):
	if places > digs or digs < 0 or places < 0:
		raise Exception('given arguments are out of definition range.')
	
	# current permutation (starting with the lowest possible combination)
	C = list(range(places))
	
	# list for all possible unique permutations
	P = [C[:]]
	
	# number of possible unique permutations without repetition (combinations)
	B = binomial(digs,places)

	# walk through all combinations
	for i in range(B-1):
		# increase the very last digit
		C[-1] += 1
		
		# and now find out, which (allowed) combination this has to be ...
		# counter has to reach places-1 to make sure all numbers fulfill
		# the conditions (the last number should automatically be correct)
		counter = 0
		# in the case anything unexpacted happens, this is the emergency exit:
		# if the algorithm exceeds 100 repetitions it stops w/o result
		save = 0
		# run through all numbers of the combination and carry over from
		# right to left if the limit of the right number is reached
		# do this, until all numbers are in range and order
		while counter < (places-1) and save < 100:
			save += 1
			counter = 0
			for n in range(places):
				# the last number (rightmost) has the biggest range,
				# the first one (leftmost) has the smallest
				if C[n] > (digs - (places - n)):
					# in case the range is exceeded, set the digit to zero
					C[n] = 0
					# and make sure all numbers are checked again
					counter = 0
				else:
					# if a successor exists ...
					if n+1 < places:
						# ... and if it is greater than the current number ...
						if C[n] < C[n+1]:
							# ... everything is fine
							counter += 1
						else:
							# if not, carray over from successor to current number
							C[n] += 1
							# and adjust successor
							C[n+1] = C[n]
							
							# the very last number does not have a successor
							if n+2 == places:
								# and thus has to be incremented here already
								C[n+1] += 1
					
		if save == 100:
			raise Exception('Sorry, your request needed to many refinement steps.')
		
		P.append(C[:])
	
	return P

def TMavailable():
	# check, if turbomole environment is set up
	try:
		cpc_path = sp.check_output('which cpc',shell=True,stderr=sp.STDOUT,universal_newlines=True)
	except sp.CalledProcessError as callerror:
		if callerror.returncode == 127:
			print('Unexpacted failure while checking the TURBOMOLE script cpc. Is the command "which" installed on your OS?')
		return False
	except FileNotFoundError:
		return False
	
	if 'TURBOMOLE' in str(cpc_path):
		return True
	
	return False
