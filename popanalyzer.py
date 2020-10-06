#! /usr/bin/python3

#################################
# popanalyzer class definition  #
#                               #
# by Fabian                     #
# 08.05.20                      #
#################################


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
class PopAnalyzerError(Exception):
	pass


class popanalyzer:
	def __init__(self,refjob,verbose=0):
		if not isinstance(refjob,jm.tmjob):
			raise PopAnalyzerError('given argument is not an instance of tmjob!')
		
		self.refjob_ = refjob
		self.popjob_ = None
		self.vrbs_lvl_ = verbose
		
		# check if there is already a population analysis in the reference job
		if len(self.refjob_.readDataGrp('$pop')) > 0:
			for f in ['dscf_pop.out','dscf.out','ridft_pop.out','ridft.out','job.last']:
				if os.path.isfile(os.path.join(self.refjob_.path_,f)):
					# if so, assign the reference job as population analysis job
					self.popjob_ = self.refjob_
					self.popout_ = os.path.join(self.refjob_.path_,f)
					break
	
	def mulliken(self,target_dir):
		if self.popjob_:
			return
		
		# define population analysis directory and check if it already exists
		pop_dir = os.path.join(self.refjob_.path_,target_dir)
		try:
			self.popjob_ = jm.tmjob(os.path.join(pop_dir,'control'))
			if self.vrbs_lvl_ > 0:
				print(' population analysis already exist in "' + str(pop_dir)  + '"')
		except jm.TMJobHandlerError:
			# check environment
			if not TMavailable():
				raise PopAnalyzerError('turbomole is not known! please set all necessary paths.')
			
			# if no already existing TM outputs can be assigned             
			# remove old pop dir (if existing)
			if os.path.isdir(pop_dir):
				shutil.rmtree(pop_dir)
			
			# copy input to new directory
			start_dir = os.getcwd()
			os.chdir(self.refjob_.path_)
			if self.vrbs_lvl_ > 1:
				print(" copying job files ...")
				os.system('cpc ' + str(pop_dir))
				print()
			else:
				os.system('cpc ' + str(pop_dir) + ' > cpc.err 2>&1')
			os.chdir(start_dir)
			
			# run population analysis
			if not os.path.isdir(pop_dir):
				raise PopAnalyzerError('error while copying turbomole files to directory "' + str(pop_dir) + '"!')
			
			
			try:
				self.popjob_ = jm.tmjob(os.path.join(pop_dir,'control'))
				self.popjob_.addToControl(['$pop'])
				self.popjob_.run(prop=True)
			except jm.TMJobHandlerError as tmerr:
				print(tmerr)
				raise PopAnalyzerError('failed to perform Mulliken population analysis!')
		
		# read population block from Mulliken analysis
		self.popout_ = self.popjob_.getOutputFile('energy')
		if not self.popout_:
			raise PopAnalyzerError('unable to read single point output!')
		
	def printSpinDensity(self):
		if not self.popjob_:
			self.mulliken('pop')
	
		print(" Mulliken Spin densities")
		print("-------------------------")
		with open(self.popout_,'r') as fh:
			reading_pop = False
			for line in fh:
				if line.strip() == '':		# skip empty lines
					continue
				
				if reading_pop:
					if '==========' in line and reading_pop:
						break
					
					print(line.strip('\n'))
						
				if 'Unpaired electrons from D(alpha)-D(beta)' in line:
					reading_pop = True

