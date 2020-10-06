#! /usr/bin/python3

###############################
# localizer class definition  #
#                             #
# by Fabian                   #
# 20.05.20                    #
###############################


# load some helpful modules
import os
import shutil
import tmjob as jm
from tools import TMavailable


# prevent stand-alone execution
if __name__ == "__main__":
	print("This class definition is not meant to be run on its own!")
	exit()


# specialized exception
class LocalizerError(Exception):
	pass


class localizer:
	def __init__(self,refjob,verbose=0):
		if not isinstance(refjob,jm.tmjob):
			raise LocalizerError('given argument is not an instance of tmjob!')
		
		self.refjob_ = refjob
		self.locjob_ = None
		self.vrbs_lvl_ = verbose
		
		# check if there are already localized orbitals in the reference job
		if len(self.refjob_.readDataGrp('$localize')) > 0:
			if os.path.isfile(os.path.join(self.refjob_.path_,'lalp')) or os.path.isfile(os.path.join(self.refjob_.path_,'lmo')):
				# if so, assign the reference job as localized job
				self.locjob_ = self.refjob_
	
	# returns a dict with keys being the LMO numbers of LMOS with predominant contribution of atoms from AtomIndices
	# and with values being dicts each assining the charge contributions among the atoms (strongly localized MOs will
	# have contributions only from one of the atoms in AtomIndices)
	def getLMOIndices(self,AtomIndices,spin=None,tol=None):
		# AtomIndices may be a list of strings or a single string;
		# each string is expected to have the form '#line' + 'element symbol', e.g. '2fe'
		if isinstance(AtomIndices,str):
			AtomIndices = [AtomIndices]
		
		# in case there are no localized orbitals yet, produce them
		if not self.locjob_:
			self.boys()		# this sets self.locjob_
		
		boys_output = self.locjob_.getOutputFile('energy')
		if not boys_output:
			raise LocalizerError('output file from orbital localization is missing in "' + str(self.locjob_.path_) + '"!')
		
		indices = {}
		with open(boys_output) as fh:
			# set defaults for search range: all localized MOs
			start_anker = 'BOYS ORBITAL LOCALISATION'
			end_anker = '=========='
			occ_tol = 0.1 if not tol else tol
			
			# adapt search range to alpha/beta shells only
			if spin == 'alpha':
				start_anker = 'ALPHA SHELLS:'
				end_anker = 'BETA SHELLS:'
			elif spin == 'beta':
				start_anker = 'BETA SHELLS:'
			
			# expected output of boys localization looks like (here refered to as InfoBlock)
			"""
			 LOCALISED MO NO. 41     diag(fock) [lmo basis] =        -0.43505
			
			   MO contributions yielding  90.00 % of density:
			      MO (col)  energy       contribution
			    44a    1   -0.40810        0.46352
			    43a    1   -0.44892        0.46241
			
			 Mulliken contributions greater than  0.1000000:
			   1fe     0.33279   0.06095   0.02034   0.25149   0.00000
			   2fe     0.33339   0.06132   0.02046   0.25160   0.00000
			   3fe     0.33355   0.06105   0.02037   0.25213   0.00000
			"""
			
			reading_LMOs = False
			tmp_num = 0
			
			for line in fh:
				if line.strip() == '':
					continue		# skip empty lines
				
				if reading_LMOs:
					if end_anker in line:
						break
				
					if 'LOCALISED MO NO.' in line:
						# InfoBlock found, set back iterator one line
						block = self.__readLMOInfoBlock(line,fh)
						try:
							tmp_num = int(block[0].split()[3])
						except ValueError:
							tmp_num = int(block[0].split()[2].split('.')[1])
						
						tmp_occ,contribs = self.__accumulateCharges(AtomIndices,block)
						
						# accept lmos with occupation close to 1.00
						if not self.locjob_.isUHF(): tmp_occ /= 2.0
						if abs(tmp_occ - 1.00) < occ_tol:
							indices[tmp_num] = contribs
				
				if start_anker in line:
					reading_LMOs = True
					fh.readline()		# skip one line (containing "...======...") to avoid premature break
		
		return indices
	
	def __readLMOInfoBlock(self,line,fh):
		block = []
		line = line.strip()
		c = 0						# emergency stop
		while not "-------" in line and c < 40:		# 40 is arbitrary and can be increased ad lib
			if line != '':
				block.append(line)
			line = fh.readline().strip()
			c += 1
		
		if len(block) == 0:
			raise LocalizerError('unable to read localization data from "' + str(fh.name) + '"!')
		
		return block
	
	def __accumulateCharges(self,idxs,block):
		contribs = {i:0.0 for i in idxs}
		occ = 0.0
		
		# find the beginning of the Mulliken analysis, add one to get in the next line
		try:
			start = block.index('Mulliken contributions greater than  0.1000000:') + 1
		except ValueError:
			raise LocalizerError('unable to locate Mulliken charge contributions in\n' + block[0])
		
		if start >= len(block):
			raise LocalizerError('Mulliken charge contribution table seems to be empty in\n' + block[0])
		
		for line in block[start:]:
			words = line.split()
			
			if words[0] in idxs:
				try:
					tmp_occ = float(words[1])
				except ValueError:
					try:
						tmp_occ = float(words[2])
					except ValueError:
						continue
				
				contribs[words[0]] = tmp_occ
				occ += tmp_occ
		
		return (occ,contribs)
	
	def findBetaElectrons(self,AtomIndices,tol=0.4):
		beta_lmos = self.getLMOIndices(AtomIndices,spin='beta',tol=tol)
		
		# pimp extracted data
	#	num_beta = len(beta_lmos)		# this statement is not true for delocalized beta e- in bonds, e.g. Fe-O
		total_contribs = {i:0.0 for i in AtomIndices}
		tmp_beta = 0.0
		for lmo,con in beta_lmos.items():
			for idx in AtomIndices:
				try:
					total_contribs[idx] += con[idx]
					tmp_beta += con[idx]
				except KeyError:
					raise LocalizerError('internal error. something went wrong while assigning the charge contributions of LMO #{}'.format(lmo))
		
		num_beta = int(round(tmp_beta))
		return (num_beta,total_contribs)
	
	def boys(self,target_dir):
		if self.vrbs_lvl_ > 1:
			print("entering BOYS LOCALIZATION")
		
		if self.locjob_:
			return
		
		# fetch necessary data
		numE = self.refjob_.getNumE()
		numVE = self.refjob_.getNumVE()
		MS = self.refjob_.getMS()
		if numE == 0 or numVE == 0:
			raise LocalizerError('problem concerning number of (valence-)electrons: #E = {:d}, #VE = {:d}'.format(numE,numVE))
			
		# open shell case (UHF)
		startMO = numE - numVE + 1
		endMO = numE
		
		# closed shell case (RHF)
		if MS == 0.0:
			startMO = numVE//2 + 1
			endMO = numE//2
		
		# define localization directory and check if it already exists
		loc_dir = os.path.join(self.refjob_.path_,target_dir)
		try:
			self.locjob_ = jm.tmjob(os.path.join(loc_dir,'control'))
			if self.vrbs_lvl_ > 0:
				print(' localized orbitals already exist in "' + str(loc_dir)  + '"')
		except jm.TMJobHandlerError:
			# check environment
			if not TMavailable():
				raise LocalizerError('turbomole is not known! please set all necessary paths.')
			
			# if no already existing TM outputs can be assigned		
			# remove old loc dir (if existing)
			if os.path.isdir(loc_dir):
				shutil.rmtree(loc_dir)
		
			# copy input to new directory
			start_dir = os.getcwd()
			os.chdir(self.refjob_.path_)
			if self.vrbs_lvl_ > 1:
				print(" copying job files ...")
				os.system('cpc ' + str(loc_dir))
				print()
			else:
				os.system('cpc ' + str(loc_dir) + ' > cpc.err 2>&1')
			os.chdir(start_dir)
		
			if self.vrbs_lvl_ > 0:
				print(" localizing the valence orbitals " + str(startMO) + "-" + str(endMO) + " ...")
				print()
		
			# run localization job
			if not os.path.isdir(loc_dir):
				raise LocalizerError('error while copying turbomole files to directory "' + str(loc_dir) + '"!')
			
			try:
				self.locjob_ = jm.tmjob(os.path.join(loc_dir,'control'))
				self.locjob_.addToControl(['$localize mo ' + str(startMO) + '-' + str(endMO)])
				self.locjob_.run(prop=True)
			except jm.TMJobHandlerError as tmerr:
				print(tmerr)
				raise LocalizerError('failed to localize orbitals!')
		
		# now, self.locjob_ should contain boys localized orbitals
		# make sure the respective files exist
		if MS == 0.0:
			lmo_file = os.path.join(self.locjob_.path_,'lmos')
		else:
			lmo_file = os.path.join(self.locjob_.path_,'lalp')
		
		if not (os.path.isfile(lmo_file) and os.stat(lmo_file).st_size > 5):
			raise LocalizerError('localization did not produce desired orbital files!')
		
