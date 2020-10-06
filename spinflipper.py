#! /usr/bin/python3

#################################
# spinflipper class definition  #
#                               #
# by Fabian                     #
# 21.05.20                      #
#################################


# load some helpful modules
import os
import shutil
import re
from itertools import combinations, product
from operator import itemgetter
import numpy as np
import tmjob as jm
from tools import TMavailable


# prevent stand-alone execution
if __name__ == "__main__":
	print("This class definition is not meant to be run on its own!")
	exit()


# Localized Molecular Orbital container class
class LMOset:
	def __init__(self,lj=None):
		self.locjob_ = lj
		self.metal_alpha_idxs_ = {}
		self.num_beta_electrons_ = 0
		self.partition_beta_electrons_ = None
		self.other_alpha_idxs_ = None


# specialized exception
class SpinFlipperError(Exception):
	pass


class spinflipper:
	def __init__(self,refjob,hs_lmos,dirprefix='',scaredy_cat=False,ox_tol=0.1,verbose=0):
		if not isinstance(refjob,jm.tmjob):
			raise SpinFlipperError('given refernce job is not an instance of tmjob!')
		
		if not isinstance(hs_lmos,LMOset):
			raise SpinFlipperError('given LMOs are not an instance of LMOset!')
		
		# check environment
		if not TMavailable():
                        raise SpinFlipperError('turbomole is not known! please set all necessary paths.')
		
		self.refjob_ = refjob
		self.lmos_ = hs_lmos
		self.prefix_ = dirprefix
		self.lsjobs_ = []
		self.vrbs_lvl_ = verbose
		
		if not self.refjob_.isUHF() or self.refjob_.getMS() < 1.0:
			raise SpinFlipperError('given input job has to be high spin!')
		
		self.ref_alpha_ = self.refjob_.getNumE('alpha')
		self.ref_beta_ = self.refjob_.getNumE('beta')
		self.ref_MS_ = self.refjob_.getMS()
		
		m = 0 if scaredy_cat else 1
		self.__fetchMOs()
		self.beta_occupations_ = self.__createBetaOccList(mode=m,ox_tol=ox_tol)
		
		# user info
		print()
		print(" Spin Flipper")
		print("--------------")
		print(" #core orbitals: {:4d}   #LMOs: {:4d} ({:4d} flipable, {:4d} others)   #virtual orbitals: {:4d}   #excess electrons: {:4d}".format(len(self.core_mos_), \
		len(self.loc_mos_), len(self.loc_mos_)-len(self.other_lmos_), len(self.other_lmos_), len(self.virt_mos_), self.lmos_.num_beta_electrons_),flush=True)
		
	def __fetchMOs(self):
		# general electronic information
		num_e = self.ref_alpha_ + self.ref_beta_
		num_ve = self.refjob_.getNumVE()
		num_core_alpha = (num_e - num_ve) // 2
		
		# gather all canonical MO coefficients
		can_mos = list(self.refjob_.getTextMOs(spin='alpha',sequential=True).values())
		self.core_mos_ = can_mos[:num_core_alpha]
		self.virt_mos_ = can_mos[self.ref_alpha_:]
		
		# get all LMO coefficients
		self.loc_mos_ = list(self.lmos_.locjob_.getTextMOs(spin='alpha',local=True,sequential=True).values())
		
		# extract LMOs of unflipable atoms
		self.other_lmos_ = [self.loc_mos_[i-1] for i in self.lmos_.other_alpha_idxs_]
		
		# remember file header of alpha and beta mo files
		with open(os.path.join(self.refjob_.path_,'alpha'),'r') as fh:
			self.ref_alpha_head_ = fh.readline()
		
		with open(os.path.join(self.refjob_.path_,'beta'),'r') as fh:
			self.ref_beta_head_ = fh.readline()
	
	# creates dicts containing meaningful integer distributions of already existing beta electrons at the flipable metal centers
	# this is necessary for reduced transition metals; example output:
	# [{'1fe':1, '3fe':0, '4co':1, '7ni':3, '9ni':2}, {'1fe':0 ,'3fe':1, '4co':1, '7ni':3, '9ni':2},
	#  {'1fe':1, '3fe':0, '4co':1, '7ni':2, '9ni':3}, {'1fe':0 ,'3fe':1, '4co':1, '7ni':2, '9ni':3}]
	# for a system containing: 1xFe(II), 1xFe(III), 1xCo(III), 1xNi(II) & 1xNi(III)
	#
	# two modi exist:
	# mode=0 is rather exhaustive, i.e. it produces all electron distributions
	# mode=1 reduces the number of distributions if the beta electrons in a domain are properly localised,
	# since they don't have to be distributed over all centers in the domain
	# e.g. Fe4O5+ --> 1xFe(II) & 3xFe(III)
	# mode=0:	[{'1fe':1, '2fe':0, '3fe':0, '4fe':0}, {'1fe':0, '2fe':1, '3fe':0, '4fe':0},
	# 		 {'1fe':0, '2fe':0, '3fe':1, '4fe':0}, {'1fe':0, '2fe':0, '3fe':0, '4fe':1}]
	# mode=1:	[{'1fe':1, '2fe':0, '3fe':0, '4fe':0}]
	def __createBetaOccList(self,mode=0,ox_tol=0.1):
		# find metal domains (i.e. kinds of metals)
		metal_centers = self.lmos_.partition_beta_electrons_.keys()					# e.g. ['1fe','3fe','4co','7ni','9ni']
		metal_domains = list(dict.fromkeys([re.findall("[a-z]+", mc)[0] for mc in metal_centers]))	# e.g. ['fe','co','ni']
		
		# count beta electrons per domain
		num_beta_per_dom = {d:0.0 for d in metal_domains}
		num_atoms_per_dom = {d:0 for d in metal_domains}
		ambig_occ_per_dom = {d:False for d in metal_domains}
		for center in metal_centers:
			for dom in metal_domains:
				if dom in center:
					occ = self.lmos_.partition_beta_electrons_[center]
					num_beta_per_dom[dom] += occ
					num_atoms_per_dom[dom] += 1
					ambig_occ_per_dom[dom] = abs(occ - round(occ)) >= ox_tol		# check for ambiguous occupations
		
		# sanity check I: compare number of beta electrons from num_beta_per_dom and self.lmos_.num_beta_electrons_
		assert int(round(np.array(list(num_beta_per_dom.values())).sum())) == self.lmos_.num_beta_electrons_, \
		"electron number not matching in {} (desired value: {})".format(num_beta_per_dom,self.lmos_.num_beta_electrons_)
		
		# produce integer occupation lists for beta electrons
		beta_distribs = []
		for dom in metal_domains:
			occ_dicts = []
			centers_in_dom = [center for center in metal_centers if dom in center]
			if mode == 0 or ambig_occ_per_dom[dom]:
				if self.vrbs_lvl_ > 1:
					print("  -> using exhaustive excess electron redistribution mode in domain {} (mode={}, ambig_occ={})".format(dom,mode,ambig_occ_per_dom[dom]))
				equal_distrib = int(round(num_beta_per_dom[dom]) // num_atoms_per_dom[dom])
				individual = int(round(num_beta_per_dom[dom]) % num_atoms_per_dom[dom])
				for comb in list(combinations(centers_in_dom,individual)):
					occ_nums = {}
					for center in centers_in_dom:
						occ_nums[center] = equal_distrib
						occ_nums[center] += 1 if center in comb else 0
					occ_dicts.append(occ_nums)
			else:	# in case the beta occupation pattern from LMOs is rather definite produce only one occ list based on those numbers
				occ_dicts.append({center:int(round(self.lmos_.partition_beta_electrons_[center])) for center in centers_in_dom})
			
			beta_distribs.append(occ_dicts)
		beta_occupations_raw = product(*beta_distribs)		# produces all combinations of dicts stored in beta_distribs (i.e. take on dict per domain)
		beta_occupations = []
		for occ in beta_occupations_raw:			# this glues together the domain-specific dicts to form the desired output
			lst = {}
			for d in occ:
				lst.update(d)
			
			# sanity check II: count electrons in each dict in beta_occupations and compare those with self.lmos_.num_beta_electrons_
			assert np.array(list(lst.values())).sum() == self.lmos_.num_beta_electrons_, \
			"electron number not matching in {} (desired value: {})".format(lst,self.lmos_.num_beta_electrons_)
			
			beta_occupations.append(lst)
		
		assert len(beta_occupations) >= 1, "not good!"
		
		return beta_occupations
	
	def __writeOrbFile(self,path,cont):
		with open(path,'w') as fh:
			for num,mo in enumerate(cont):
				if num == 0:			# write header
					fh.write(mo)
					continue
				
				orb_num = int(mo[0].split()[0])
				
				if not orb_num == (num):
					mo[0] = mo[0].replace(str(orb_num),str(num),1)
				
				for line in mo:			# write MO coeffs
					fh.write(line)
			fh.write('$end')
	
	def __createLSJob(self,dir_name):
		flip_dir = os.path.join(self.refjob_.path_,dir_name)
		
		if os.path.isdir(flip_dir):
			shutil.rmtree(flip_dir)
		
		# copy input to new file
		start_dir = os.getcwd()
		os.chdir(self.refjob_.path_)
		if self.vrbs_lvl_ > 1:
			print("copying job files ...")
			os.system('cpc ' + str(dir_name))
			print()
		else:
			os.system('cpc ' + str(dir_name) + ' > cpc.err 2>&1')
		
		os.chdir(start_dir)
		
		if not os.path.isdir(flip_dir):
			raise SpinFlipperError('error while copying turbomole files to directory "' + str(flip_dir) + '"!')
		
		# set up new TM job
		try:
			lsjob = jm.tmjob(os.path.join(flip_dir,'control'))
		except jm.TMJobHandlerError as tmerr:
			print(tmerr)
			raise SpinFlipperError('unable to set up low spin job!')
		
		alpha_path = os.path.join(lsjob.path_,'alpha')
		beta_path = os.path.join(lsjob.path_,'beta')
		
		os.remove(alpha_path)
		os.remove(beta_path)
		
		return lsjob
	
	def flip(self,centers):
		if self.vrbs_lvl_ > 1:
			print("entering FLIPPING")
		
		# check input
		if len(centers) == 0:
			raise SpinFlipperError('no centers specified to be flipped!')
		
		if self.vrbs_lvl_ >= 0:
			print(" Flipping electrons at " + ", ".join(centers))
		if self.vrbs_lvl_ > 0 and len(self.beta_occupations_) > 1:
			print("  There are {:d} possibilities to distribute the existing excess electrons among the metal centers.".format(len(self.beta_occupations_)))
		
		num_existing_lsjobs = len(self.lsjobs_)
		
		# calculate sorting weights of LMOS depending on occupation
		for nr,beta_occ in enumerate(self.beta_occupations_):
			lmo_weights_alpha = []
			lmo_weights_beta = []
			num_flip_alpha = 0
			for center,lmo_idxs in self.lmos_.metal_alpha_idxs_.items():
				for c,i in enumerate(lmo_idxs):
					if center in centers:
						w_a = 10 if c < beta_occ[center] else 100
						w_b = 0
						num_flip_alpha += 1
					else:
						w_a = 0
						w_b = 10 if c < beta_occ[center] else 100
					lmo_weights_alpha.append((i,w_a))
					lmo_weights_beta.append((i,w_b))
			
			# sort LMO indices of metal centers
			alpha_lmo_idxs = sorted(lmo_weights_alpha,key=itemgetter(1))
			beta_lmo_idxs = sorted(lmo_weights_beta,key=itemgetter(1))
			
			alpha_metal_lmos = [self.loc_mos_[i-1] for i,w in alpha_lmo_idxs]
			beta_metal_lmos = [self.loc_mos_[i-1] for i,w in beta_lmo_idxs]
			
			# glue MO coeffs together in the right way
			new_alpha_mos = [self.ref_alpha_head_] + self.core_mos_ + self.other_lmos_ + alpha_metal_lmos + self.virt_mos_
			new_beta_mos = [self.ref_beta_head_] + self.core_mos_ + self.other_lmos_ + beta_metal_lmos + self.virt_mos_
			
			# calculate new alpha and beta occupation numbers
			add_alpha = [iw[1] for iw in lmo_weights_alpha].count(10)		# count the amount of weight 10
			add_beta = [iw[1] for iw in lmo_weights_beta].count(10)			# which stands for excess beta electrons
			assert add_alpha+add_beta == self.lmos_.num_beta_electrons_, \
			"incorrect number of excess electrons: {} a, {} b; (desired value: {})".format(add_alpha,add_beta,self.lmos_.num_beta_electrons_)
			new_alpha_occ = self.ref_alpha_ - num_flip_alpha + add_alpha
			new_beta_occ  = self.ref_beta_  + num_flip_alpha - add_alpha
			new_MS = abs(new_alpha_occ - new_beta_occ) / 2.0
			
			# inform the user
			if self.vrbs_lvl_ > 0:
				print(" {:4d}.".format(nr+1), end='')
				for center in beta_occ:
					print(" {}: {:3d};".format(center, beta_occ[center]), end='')
				print(flush=True)
				print("       flipped alpha electrons: {}   flipped excess electrons: {}".format(num_flip_alpha, add_alpha))
			
			# creat new directory
			# the directory's name consists of the label of the flipped centers ...
			dir_name = self.prefix_
			for item in centers:
				dir_name += str(item)
			# ... the new spin multiplicity ...
			dir_name += '_' + str(int(2*new_MS+1)) + 'tet'
			# ... and a running number
			if len(self.beta_occupations_) > 1: dir_name += '_' + str(nr+1)
			
			self.lsjobs_.append(self.__createLSJob(dir_name))
			
			# create new orbital files
			if self.vrbs_lvl_ > 1:
				print('creating new orbital files ...')
			
			# write out new alpha and beta orbitals
			self.__writeOrbFile(os.path.join(self.lsjobs_[-1].path_,'alpha'),new_alpha_mos)
			self.__writeOrbFile(os.path.join(self.lsjobs_[-1].path_,'beta'),new_beta_mos)
			
			# update control file
			try:
				iter_limit = int(self.lsjobs_[-1].readDataGrp('$scfiterlimit')[0].split()[1])
			except (IndexError,ValueError):
				iter_limit = 0
			new_iter_limit = iter_limit if iter_limit > 500 else 500
			self.lsjobs_[-1].removeFromControl(['$alpha shells','$beta shells','$scfdamp','$scforbitalshift','$scfiterlimit'])
			self.lsjobs_[-1].addToControl([	'$alpha shells',' a       1-{:d}                     (1)'.format(new_alpha_occ),\
							'$beta shells', ' a       1-{:d}                     (1)'.format(new_beta_occ),\
							'$scfdamp   start=5.500  step=0.050  min=0.500',\
							'$scforbitalshift  automatic=1.0',\
							'$scfiterlimit {:d}'.format(new_iter_limit)])
			self.lsjobs_[-1].updateControl()
			
			assert self.lsjobs_[-1].getMS() < self.ref_MS_, \
			"new low spin job {} has wrong occupation, old MS: {}, new MS: {}".format(self.lsjobs_[-1],self.ref_MS_,self.lsjobs_[-1].getMS())
			
			# run job!
			# Attention! This is only needed for testing reasons!
			#try:
			#	self.lsjobs_[-1].run()
			#except jm.TMJobHandlerError as tmerr:
			#	print(tmerr)
			#	raise SpinFlipperError('failed to flip spins! please check job at {}'.format(self.lsjobs_[-1]))
		
		assert len(self.lsjobs_) == num_existing_lsjobs + len(self.beta_occupations_), \
		"there is not the correct number of low spin jobs, old: {:d}, new {:d}, supposed added: {:d}".format(num_existing_lsjobs, len(self.lsjobs_), len(self.beta_occupations_))
		
		return self.lsjobs_[-len(self.beta_occupations_):]

