#! /usr/bin/python3

###############################
# lowSpin  v1.3.0             #
# (23.05.20)                  #
#                             #
# by Fabian                   #
# 30.11.15                    #
###############################


# check python version
import sys
if not sys.version_info.major == 3 and sys.version_info.minor >= 4:
	print("This program requires at least python 3.4!")
	exit()

# load some helpful modules
import argparse
import os
import glob
import tmjob as jm
import QueueSys as qs
import localizer as lc
import popanalyzer as pa
import spinflipper as sf
import PSE
import re
import math
from itertools import combinations
#from tools import getCombinations
#import traceback
from datetime import datetime


# global definition
ox_state = {-5:'-V', -4:'-IV', -3:'-III', -2:'-II', -1:'-I', 0:'0', 1:'I', 2:'II', 3:'III', 4:'IV', 5:'V', 6:'VI', 7:'VII', 8:'VIII'}

#############################################
# Spin flipping algorithm                   #
#############################################

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Semi-automatic spin flipping algorithm for TURBOMOLE')
	parser.add_argument('hsjob',metavar='JOB',nargs='?',default='.',help='relative path to the high spin job (default: cwd)')
	parser.add_argument('--pbs-script','-p',nargs=1,metavar='SCRIPT',help='relative path to a PBS job script')
	parser.add_argument('--scaredy-cat','-s',action='store_true',help='use the exhaustive algorithm when redistributing the excess electrons and produce a lot of low spin systems')
	parser.add_argument('--analysis','-a',action='store_true',help='performe only localization and analysis')
	parser.add_argument('--beta-tolerance','-b',nargs=1,metavar='TOL',type=float,default=[0.4],help='Accepted occupation deviation from 1.0 for beta LMOs when searching excess electrons (default: 0.4)')
	parser.add_argument('--alpha-tolerance','-t',nargs=1,metavar='TOL',type=float,default=[0.1],help='Accepted occupation deviation from 1.0 for alpha LMOs when searching flipable electrons (default: 0.1)')
	parser.add_argument('--ox-tolerance','-o',nargs=1,metavar='TOL',type=float,default=[0.1],help='Accepted deviation from proper integer occupation number for determination of the oxidation state (default: 0.1)')
	parser.add_argument('--verbose','-v',nargs=1,metavar='LEVEL',type=int,default=[0],help='change verbose level (0 means off)')
	
	args = parser.parse_args()
	
	hs_job_path = os.getcwd()
	if not args.hsjob is None:
		if os.path.isdir(os.path.abspath(args.hsjob)):
			hs_job_path = os.path.abspath(args.hsjob)
	
	if not args.analysis:
		try:	
			pbs_script_path = os.path.abspath(args.pbs_script[0]) if args.pbs_script else glob.glob("*.job")[0]
		except IndexError:
			print("Please provide a PBS job script either via the -p option or as .job file in the cwd.")
			exit()
	
	vrbs_level = args.verbose[0]
	
#	print("0: ", datetime.now().time(),flush=True)
	print()
	print('\t+----------------------+')
	print('\t|    lowSpin v1.3.0    |')
	print('\t+----------------------+')
	print(flush=True)
		
	try:
		# set up environment
		highspinjob = jm.tmjob(os.path.join(hs_job_path,'control'))
		
		# analyze input job
		element_abundance = highspinjob.getElementAbundances()	
	
		# possible candidates for spin flipping ... this list can be expanded on demand
		candidates = ['cr','mn','fe','co','ni','cu']
		
		# create a list of all present metal atoms (those that can be flipped) in the reference job
		metal_centers = []
		for metal in candidates:
			metal_centers += highspinjob.getAtomIndexList(metal)	# appends an empty list, if the resp. metal is not in the molecule
										# metal_centers is now a list like ['9cr','2fe','3fe','5co']
		
		# get number of metal centers
		num_centers = len(metal_centers)
		if num_centers < 2 and not args.analysis:
			print("there is only one metal atom in the system; this program can't help you here.")
			exit()
		
		# fetch general electronic information and spin density population analysis
		popanalyzer = pa.popanalyzer(highspinjob,vrbs_level)
		popanalyzer.mulliken('pop')
		
		# localize MOs (if not already done)
		localizer = lc.localizer(highspinjob,vrbs_level)
		localizer.boys('loc')
		
		# collect all necessary data from the (localized) high spin system
		hs_lmos = sf.LMOset(localizer.locjob_)		# container for lmo infos needed for spin flipping
		print(" Evaluating Localized MOs ...\n  ",end="",flush=True)
		for center in metal_centers:
			# dict like {'1fe':[1,2,4,5],'2fe':[3,7,8,9,10],'4co':[23,24,26,30]}
			hs_lmos.metal_alpha_idxs_[center] = localizer.getLMOIndices(center,spin='alpha',tol=args.alpha_tolerance[0]).keys()
			print("#",flush=True,end="")
		print("\n",flush=True)
		# for spin flipping only the number of already existing beta electrons in the valence states of the metals
		# that will be flipped are of interest, not their actual LMOs since they are assumed to be delocalized anyways
		hs_lmos.num_beta_electrons_,hs_lmos.partition_beta_electrons_ = localizer.findBetaElectrons(metal_centers,tol=args.beta_tolerance[0])
		# for reordering the LMOs it is necessary to know the remainder of orbitals that are not located at one of the
		# metal sites but belong to the valence region
		num_alpha_VE = highspinjob.getNumE('alpha') - (highspinjob.getNumE() - highspinjob.getNumVE()) // 2
		metals_idxs = [idx for i in hs_lmos.metal_alpha_idxs_.values() for idx in i]
		hs_lmos.other_alpha_idxs_ = [idx for idx in range(1,num_alpha_VE+1) if not idx in metals_idxs]
		
		# give a summary of the high spin reference state
		popanalyzer.printSpinDensity()
		
		print(flush=True)
		print()
		print(" Analysis of Localized MOs")
		print("---------------------------")
		print()
		print("  center   #alpha   #beta   ox. state")
		print(" -------------------------------------")
		for center in hs_lmos.metal_alpha_idxs_:
			a = len(hs_lmos.metal_alpha_idxs_[center])
			b = hs_lmos.partition_beta_electrons_[center]
			e = re.findall("[a-z]+", center)[0]
			try:
				d = PSE.VE[e] - a - b
			except IndexError:
				print("  The element {} is not listed in the program's periodic table so far... :(".format(e))
				ox = "---"
			else:
				try:
					if abs(d - round(d)) < args.ox_tolerance[0]:
						ox = ox_state[int(round(d))]
					else:
						ox = ox_state[math.floor(d)] + "/" + ox_state[math.ceil(d)]
				except IndexError:
					print("  Oxidation state {:f} of element {} seems strange.".format(d,e))
					ox = "---"
			
			print("  {:^6}     {:3d}     {:4.2f}   {:^9}".format(center,a,b,ox))
		print(" -------------------------------------")
		print(flush=True)
		
		if not args.analysis:
			# set up queuing system
			submitter = qs.QueueSys(pbs_script_path)
			
			# set up the spin flipper (tell it to create new subdirs beginning with the 'flip_' and use the exhaustive algorithm on user request)
			spinflipper = sf.spinflipper(highspinjob,hs_lmos,'flip_',args.scaredy_cat,args.ox_tolerance[0],vrbs_level)
			
			# loop through all combinations
			combinations = combinations(list(range(num_centers)),int(num_centers//2))
			for comb in combinations:
				# associate the metal centers to the numbers in the current combination
				flip_centers = [metal_centers[item] for item in comb]
				
				# produce new job(s) with the current centers spin flipped
				lowspinjobs = spinflipper.flip(flip_centers)
				
				# ask user whether really to start the job
				q_start = input("  -> Submit this batch of {:d} job(s)? (default: yes)> ".format(len(lowspinjobs))).lower() in ['n','no','0']
				if q_start:
					q_keep = input("  -->  Keep job files? (default: no)> ").lower() in ['y','yes','1']
				
				for ls in lowspinjobs:
					if q_start:
						if not q_keep:
							ls.remove()
					else:
						submitter.schedule(ls.path_)
	
	except jm.TMJobHandlerError as tmerr:
		print("Error while evaluating TM job data:")
		print(tmerr)
		#traceback.print_exc()
		exit()
	except lc.LocalizerError as locerr:
		print("Error while localizing orbitals:")
		print(locerr)
		#traceback.print_exc()
		exit()
	except pa.PopAnalyzerError as paerr:
		print("Error while evaluating orbital population:")
		print(paerr)
		#traceback.print_exc()
		exit()
	except qs.QueueSysError as qserr:
		print("Error while operating with the queuing system:")
		print(qserr)
		#traceback.print_exc()
		exit()
	except sf.SpinFlipperError as sferr:
		print("Error while creating low spin job:")
		print(sferr)
		#traceback.print_exc()
		exit()
	except:
		print("Unexpected error:")
		raise
		exit()
	
	print()
	print("Done!")

