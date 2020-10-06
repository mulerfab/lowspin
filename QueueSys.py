#! /usr/bin/python3

###############################
# modified copy               #
# python3 version             #
###############################

###############################
# multiPot  v1.0.0            #
#                             #
# by Fabian                   #
# 02.01.18                    #
#                             #
# module definition for:      #
# QueueSys                    #
###############################


import os
from shutil import copyfile
import subprocess as sp


# prevent stand-alone execution
if __name__ == "__main__":
	print("This class definition is not meant to be run on its own!")
	exit()


# specialized exception
class QueueSysError(Exception):
	pass


# class definition
class QueueSys:
	def __init__(self, job_script):
		if not os.path.isfile(job_script):
			raise QueueSysError('The job script file "' + str(job_script) + '" does not exist!')
		
		self.job_script_ = os.path.abspath(job_script)
		
		self.stat_cmd_ = 'qstat'		# local:	['locque.sh', '-l']
		self.sub_cmd_ = 'qsub'			#		['locque.sh']
		
		try:
			output = sp.check_output(self.stat_cmd_,shell=True,stderr=sp.STDOUT,universal_newlines=True)
		except sp.CalledProcessError:
			raise QueueSysError('No queuing system found!')
		
		self.user_ = sp.check_output('whoami',shell=True,stderr=sp.STDOUT,universal_newlines=True).strip()
	
	
	def get_schedule(self):
		# fetch the status table from the queuing system as continuous string
		raw_stat = sp.check_output(self.stat_cmd_,shell=True,stderr=sp.STDOUT,universal_newlines=True)
		
		if raw_stat.strip() == '':
			return []
		
		# split the string in lines and the lines in words, discard first two lines (table head)
		# the order per line is then: job ID, job name, user, used cpu time, status, queue
		stat = [line.split() for line in raw_stat.strip().split('\n')[2:]]
		
		# sanity check
		if len(stat) != 0:
			if len(stat[0]) != 6:
				raise QueueSysError('Unable to read job information from command ' + str(self.stat_cmd_))
		
		# post processing...
		sched = []
		for line in stat:
			if line[2] != self.user_:		# remove all jobs from other users
				continue
			
			if line[4].strip() not in 'QRECH':	# remove unknown statuses
				print('  -> Unexpacted status of job ' + str(line[0]) + ' with name ' + str(line[1]))
				continue
			
			line[0] = int(line[0].split('.')[0])	# remove the cluster name from the job ID
			
			sched.append(line)
		
		return sched
	
	
	def __status_validity_check(self, stat, allowed, funcName):
		try:
			l = len(stat)
		except TypeError:
			raise QueueSysError('Unexpacted status in ' + str(funcName) + ': argument needs to be a string!')
		
		if l > (len(allowed) - 1):
			raise QueueSysError('Unexpacted status in ' + str(funcName) + ': argument must not be longer than ' + str(len(allowed)-1) + ' letters!')
		
		status = list(stat)
		for s in status:
			if s not in allowed:
				raise QueueSysError('Unexpacted status in ' + str(funcName) + ': ' + s + ' is not known!')
	
	
	def get_num_jobs(self, status):
		self.__status_validity_check(status, 'QRECH', 'get_num_jobs')
		
		counter = 0
		for line in self.get_schedule():
			if line[4] in status:
				counter += 1
		
		return counter
	
	
	def get_job_IDs(self, status):
		self.__status_validity_check(status, 'QRECH', 'get_job_IDs')
		
		ids = []
		for line in self.get_schedule():
			if line[4] in status:
				ids.append(line[0])
		
		return ids
	
	
	def schedule(self, job_path):
		if not os.path.isdir(job_path):
			raise QueueSysError('Unable to submit job from "' + str(job_path) + '". Path does not exist.')
		
		# copy the given pbs script file to the directory of the new job
		script_name = os.path.basename(self.job_script_)
		new_script_path = os.path.join(job_path, script_name)
		copyfile(self.job_script_, new_script_path)
		
		# adopt the top most directory name as job name
		job_name = os.path.basename(os.path.normpath(job_path))
		
		# insert the PBS line for the job name in the job script
		script_cont = []
		with open(new_script_path, "r") as fh:
			for line in fh:				# lines are with trailing new-line
				if '#PBS -N' in line:
					continue
				
				script_cont.append(line)
				
				if '#!' in line:
					script_cont.append('#PBS -N ' + str(job_name) + '\n')
		
		with open(new_script_path, "w") as fh:
			for line in script_cont:
				fh.write(line)
		
		# submit the job to the queuing system
		cwd = os.getcwd()
		os.chdir(job_path)
		try:
			job_id = sp.check_output(self.sub_cmd_ + " " + script_name,shell=True,stderr=sp.STDOUT,universal_newlines=True)
		except sp.CalledProcessError:
			raise QueueSysError('Unable to submit job "' + str(new_script_path) + '". The command ' + str(self.sub_cmd_) + ' failed.')
		os.chdir(cwd)
		
		job_id = int(job_id.split('.')[0])
		
		return job_id
	

