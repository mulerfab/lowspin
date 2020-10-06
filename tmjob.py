#! /usr/bin/python3

###############################
# tmjob class definition      #
#                             #
# by Fabian                   #
# 21.05.20                    #
###############################


# load some helpful modules
import os
import shutil
import subprocess as sp
import re
from copy import deepcopy
import PSE


# prevent stand-alone execution
if __name__ == "__main__":
	print("This class definition is not meant to be run on its own!")
	exit()


# specialized exception
class TMJobHandlerError(Exception):
	pass


#############################################
# Turbomole Job Handler                     #
#############################################

class tmjob:
	def __init__(self,control_path):
		if not os.path.isfile(control_path):
			raise TMJobHandlerError('control file at "' + str(control_path) + '" does not exist!')
		
		self.control_path_ = control_path
		self.path_ = os.path.dirname(control_path)
		
		self.element_abundances_ = None
		self.num_e_ = None
		
		with open(self.control_path_,'r') as fh:
			self.control_ = fh.readlines(1024*1024)
	
	def __repr__(self):
		return str(self.control_path_)
	
	def __str__(self):
		return str(self.path_)
	
	def updateControl(self):
		os.remove(self.control_path_)
		
		with open(self.control_path_,'w') as fh:
			for line in self.control_:
				if line.strip() != '':
					fh.write(line.strip('\n') + '\n')
	
	def addToControl(self,cmd_lines):
		if not isinstance(cmd_lines,list):
			cmd_lines = [cmd_lines]
		
		new_control = []
		for line in self.control_:
			if '$end' in line:
				for cmd in cmd_lines:
					new_control.append(cmd)
			
			new_control.append(line)
		
		self.control_ = new_control
	
	def replaceControlLine(self,old_line,new_line):
		new_control = []
		for line in self.control_:
			if old_line in line:
				new_control.append(new_line)
			else:
				new_control.append(line)
		
		self.control_ = new_control
	
	def removeFromControl(self,cmd_lines):
		new_control = []
		
		remove = False
		for line in self.control_:
			if '$' in line and remove:
				remove = False
			
			for cmd in cmd_lines:
				if cmd in line:
					remove = True
			
			if not remove:
				new_control.append(line)
		
		self.control_ = new_control
	
	def readDataGrp(self,grp_key,ext_file=None,strip=True):
		data = self.control_
		if ext_file:
			with open(os.path.join(self.path_,ext_file),"r") as fh:
				data = fh.readlines(1024*1024*500)
		
		reading_data = False
		grp = []
		for line in data:
			if '#' in line.split()[0]:
				continue
			
			if reading_data:
				if '$' in line:
					break
				else:
					grp.append(line.strip() if strip else line)
			
			if grp_key in line:
				if len(line.split()) > 1 and 'file' in line.split()[1] and not ext_file:	# to avoid infinit recursions
					return self.readDataGrp(grp_key,ext_file=line.split('=')[1].strip(),strip=strip)
				else:
					grp.append(line.strip() if strip else line)
					reading_data = True
		
		return grp
	
	def isUHF(self):
		if len(self.readDataGrp('$uhf')) == 0:
			return False
		
		return True
	
	def isC1(self):
		symm_grp = self.readDataGrp('$symmetry')
		if symm_grp[0].split()[1].lower() == 'c1':
			return True
		
		return False
	
	def isRI(self):
		if len(self.readDataGrp('$ri')) == 0:
			return False
		
		return True
	
	def getElementAbundances(self):
		if self.element_abundances_:
			return self.element_abundances_
		
		coords = self.readDataGrp('$coord')
		element_abundances = {}
		for line in coords[1:]:
			elem = line.split()[3].lower()
			if elem in element_abundances:
				element_abundances[elem] += 1
			else:
				element_abundances[elem] = 1
		
		self.element_abundances_ = element_abundances
		return element_abundances
	
	def getNumE(self,spin=None):
		if not spin and self.num_e_:
			return self.num_e_
		
		num_e = 0
		
		if spin == 'alpha' or spin == 'beta' or spin == 'closed':
			# get the whole occupation block
			occ_grp = self.readDataGrp('$' + str(spin) + ' shells')
			
			if len(occ_grp) < 2:
				raise TMJobHandlerError('unable to read orbital occupations!')
			
			# walk through the block line by line
			# a line can look like:
			#  b1      1- 12, 14                               ( 1 )
			for line in occ_grp[1:]:
				# remove the irrep and the actual occ number
				try:
					start = line.strip().index(" ")
					end = line.strip().index("(")
				except ValueError:
					raise TMJobHandlerError('format of orbital occupations unknown!')
				
				if end - start <= 0: raise TMJobHandlerError('strange formating of orbital occupations!')
				
				raw_numbers = line.strip()[start:end].strip()
				
				# split up ranges (separated by ,)
				ranges = raw_numbers.split(",")
				for r in ranges:
					# split up numbers in ranges (separated by -)
					numbers = r.strip().split("-")
					if len(numbers) != 0 and numbers[0].strip() == '': raise TMJobHandlerError('empty orbital occupation:\n{}'.format(line))
					if len(numbers) > 1:
						try:
							num_e += int(numbers[1].strip())-(int(numbers[0].strip())-1)
						except ValueError:
							raise TMJobHandlerError('integer number expected, found {}'.format(numbers))
					# or if there is only one number count it as one e-
					else:
						num_e += 1
			
			if spin == 'closed': num_e *= 2
			
	#		if not self.isC1():
	#			raise TMJobHandlerError('electron counting implemented only in C1 symmetry!')
	#		
	#		# split with regular expression:
	#		# \s means any white space, * means 0 or more times, + means 1 or more times, | means or, - means -
	#		numbers = re.split('\s*-\s*|\s+',occ_grp[1].strip())
	#		# " a   1- 57              ( 1 )" will be splitted into
	#		# ['a', '1', '57', '(', '1', ')']
	#		
	#		if len(numbers) > 2:
	#			try:
	#				num_e = int(numbers[2])
	#			except ValueError:
	#				num_e = 1
	#		else:
	#			raise TMJobHandlerError('unable to extract orbital occupation from control line!')
		else:
			if self.isUHF():
				num_alpha = self.getNumE('alpha')
				num_beta = self.getNumE('beta')
				num_e = num_alpha + num_beta
				self.num_e_ = num_e
			else:
				num_e = self.getNumE('closed')
		
		return num_e
	
	def getNumVE(self):
		numVE = 0
		for elem,abund in self.getElementAbundances().items():
			try:
				VE = PSE.VE[elem] * abund
			except KeyError:
				raise TMJobHandlerError('the element "' + str(elem) + '" is not listed in the PSE so far. Sorry.')
			
			numVE += VE
		
		charge = self.getCharge()
		return (numVE - charge)
	
	def getCharge(self):
		charge = 0
		
		charge_grp = self.readDataGrp('$charge')
		try:
			charge = int(round(float(charge_grp[1].split()[0])))
		except ValueError:
			pass
		except IndexError:
			raise TMJobHandlerError('unable to read charge from control file!')
		
		return charge
	
	def getMS(self):
		ms = 0.0
		
		if self.isUHF():
			num_alpha = self.getNumE('alpha')
			num_beta = self.getNumE('beta')
		
			ms = abs(num_alpha - num_beta) / 2.0
		
		return ms
	
	def getAtomIndexList(self,element):
		# returns a list containing an entry for each entity of "element" in the coord-file
		# the entries have the form '#line' + 'element', e.g. '3fe' for the iron atom from
		# the thrid line
		atom_index_list = []
		coords = self.readDataGrp('$coord')
		
		if len(coords) < 2:
			raise TMJobHandlerError('unable to read coordinates!')
		
		for i,line in enumerate(coords[1:]):
			if line.split()[3].lower() == element.lower():
				atom_index_list.append(str(i+1) + line.split()[3])
		
		return atom_index_list
	
	def getTextMOs(self,spin='alpha',local=False,sequential=False):
		if self.isUHF():
			if spin == 'beta':
				mos_raw = self.readDataGrp('$lmo_beta','lbet',False) if local else self.readDataGrp('$uhfmo_beta',strip=False)
			else:
				mos_raw = self.readDataGrp('$lmo_alpha','lalp',False) if local else self.readDataGrp('$uhfmo_alpha',strip=False)
		else:
			mos_raw = self.readDataGrp('$lmo','lmos',False) if local else self.readDataGrp('$scfmo',strip=False)	# NOT TESTED YET!
		
		if len(mos_raw) < 2:
			raise TMJobHandlerError('unable to read MOs!')
		
		mos = {}
		tmp_mo = []
		mo_label = 0 if sequential else ''
		for line in mos_raw[1:]:
			# search mo files for lines like:
			# "     1  a      eigenvalue=-.25612783293457D+03   nsaos=434"
			if 'eigenvalue=' in line:
				if tmp_mo:
					mos[mo_label] = deepcopy(tmp_mo)	# mo_label always holds the last found label
					tmp_mo = []
				
				if sequential:			# sequential numbering
					mo_label += 1
				else:				# labeling by irrep
					mo_label = ''.join(line.split()[:2])
					
			tmp_mo.append(line)
		
		# append very last MO
		mos[mo_label] = deepcopy(tmp_mo)
		
		return mos
	
	def getFloatMOs(self,local=False,sequential=False):
		# TBD
		pass
	
	def run(self,prop=False,opt=False,opt_flags=[],freq=False):
		self.updateControl()
		start_dir = os.getcwd()
		
		# change to job directory
		if not start_dir == self.path_:
			os.chdir(self.path_)
		
		# run single point in any case
		energy_in = 'ridft > ridft.out' if self.isRI() else 'dscf > dscf.out'
		
		# if desired calculate properties only
		if prop:
			energy_in = 'ridft -proper > ridft.out' if self.isRI() else 'dscf -proper > dscf.out'
	
		try:
			# run dscf/ridft in any case
			energy_out = sp.check_output(energy_in,shell=True,stderr=sp.STDOUT,universal_newlines=True)
			if 'abnormally' in str(energy_out):
				raise TMJobHandlerError('error while executing single point calculation in "' + str(self.path_) + '"!')
			
			# run optimization if requested
			if opt:
				jobex_flags = opt_flags
				if self.isRI() and not '-ri' in jobex_flags:
					jobex_flags += ['-ri']
				jobex_in = 'jobex ' + ' '.join(jobex_flags) + ' > jobex.out'
				jobex_out = sp.check_output(jobex_in,shell=True,stderr=sp.STDOUT,universal_newlines=True)
				if 'abnormally' in str(jobex_out):
					raise TMJobHandlerError('error while executing jobex in "' + str(self.path_) + '"!')
			
			# run frequency calculation if requested
			if freq:
				force_out = sp.check_output('aoforce > aoforce.out',shell=True,stderr=sp.STDOUT,universal_newlines=True)
				if 'abnormally' in str(force_out):
					raise TMJobHandlerError('error while executing aoforce in "' + str(self.path_) + '"!')
		except sp.CalledProcessError as tmerr:
			raise TMJobHandlerError("error while running TURBOMOLE:\nreturn code was {}\ncommand was {}".format(tmerr.returncode,tmerr.cmd))
		
		os.chdir(start_dir)
		return True
	
	def getOutputFile(self,jobtype):
		if jobtype == 'energy':
			for f in ['dscf.out','ridft.out','job.last']:
				tmp_path = os.path.join(self.path_,f)
				if os.path.isfile(tmp_path):
					return os.path.abspath(tmp_path)
		
		if jobtype == 'opt':
			tmp_path = os.path.join(self.path_,'job.last')
			if os.path.isfile(tmp_path):
				return os.path.abspath(tmp_path)
		
		if jobtype == 'freq':
			tmp_path = os.path.join(self.path_,'aoforce.out')
			if os.path.isfile(tmp_path):
				return os.path.abspath(tmp_path)
		
		return None
	
	def remove(self):
		shutil.rmtree(self.path_)
		return True
	
