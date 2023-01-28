'''
=head1 NAME

STORIcontrol.py

=head1 SYNOPSIS

	python ./STORIcontrol.py
	STORI>start <run-name> <scratch/dir> <taxa file> <windowSize> <finalMaxFams>
	STORI>start hemoglobin_eumetazoa_1x_STORI /home/ec2-user/STORI/Runfiles eumetazoa 4 20
	STORI>start hemoglobin_euk_8x_STORI /home/ec2-user/STORI/Runfiles eukaryota 4 20
	STORI>stop hemoglobin_euk_8x_STORI
	STORI>exit

=head1 DESCRIPTION

This script controls the rest of STORI, except for
STORIstats, which is independent of STORIcontrol.

STORIcontrol allows the user to start and stop protein family retrievals.
Each retrieval ("run") consists of two Perl processes, ie starting a new
run will create 2 new PIDs.

Most of the computation occurs in STORI.pl, because this is the script that
executes BLASTP searches repeatedly. Each run comprises two independent
instances of STORI.pl, one in the first PID, and the other in the second PID.

Here is an
overview of the different scripts, to help explain the purpose of STORIcontrol.py:

=head2 STORIcontrol.py accepts user input repeatedly

=over 3

=item * If the user commands "start <arguments>"
	
	Then STORIcontrol.py saves a record of the new run, and begins the run
	by invoking beginSTORI.pl
	
	beginSTORI.pl guides the user through choosing seed proteins
	and passes two randomized sets of this information to two new
	STORI.pl instances that it kicks off.
			
=item * If the user commands "stop <run name>"
	
	Then STORIcontrol.py kills the PIDs for this run.
	NB funny things may occur if the user stops runs that failed
	to start.
		
=item * If the user does nothing
	
	Then after 15 seconds STORIcontrol will cycle through its record of
	existing runs, and update the status of each run using checkSTORI.pl.
	The run parameters are saved in job_data_STORI.txt
	
	If both PIDs in a run have completed
	
		Then STORIcontrol will examine the last three agreement scores for the
		run to determine whether or not the run has converged. (Each PID could
		be considered a Markov chain, and when the chains are sufficiently similar,
		then the run has converged.)
		
		If the run converged, then it is labelled "converged" and nothing more happens
		with it.
		
		If the run has not converged, then STORIcontrol will create two new instantiations of
		STORI.pl as above, using the intersecting result of the jobs that just finished. 
		A pair of STORI.pl processes has a wall clock limit specified by the global variable
		maxRuntime, below. If the STORI.pl processes do not end of their own accord, and they
		have both exceeded maxRuntime, STORIcontrol will wait until each instance has output
		at least one GI table, then kill these processes and check for convergence.
		
		
	Otherwise
	
		STORIcontrol simply updates the agreement score between the sequence groupings,
		which is an output of checkSTORI.pl. If checkSTORI is checking a run that began
		recently, it is possible that one or both of the runs have not yet generated any
		output GI table, in which case the agreement score = -1. STORIstats will only be able
		to summarize runs whose agreement score is > 0.

=back

=cut

my $help;
GetOptions(
	   'h'             => \$help,
	   );

if( $help ) {
    exec('perldoc',$0);
    exit(0);
}
'''

import subprocess
import re
from collections import defaultdict
import time
import threading
import Queue
import os.path

home = "/home/josh"
arrayFile = home + "/STORI/job_data_STORI.txt"
checkSTORIpath = home + "/STORI/checkSTORI.pl"
continueSTORIpath_t = home + "/STORI/continueSTORIfast_t.pl"
beginSTORIpath_t = home + "/STORI/beginSTORI.pl"
maxRuntime = 2  	#maximum wall time of a STORI.pl process in hours

subprocess.check_call(["touch", arrayFile])
arrayFileBackup = home + "/STORI/job_data_STORI_backup.txt"
subprocess.check_call(["cp", arrayFile, arrayFileBackup])



class Run:
	def __init__(self, runName, parentDir, taxaFile, windowSize, finalMaxFams, STORIcount):
		self.runA_id=-1
		self.runB_id=-1
		self.runName = runName
		self.parentDir = parentDir
		self.taxaFile = taxaFile
		self.windowSize = int(windowSize)
		self.finalMaxFams = int(finalMaxFams)
		self.STORIcount = int(STORIcount)
		self.converged = False
		self.paused = False
		self.runStats = {}
		self.runStats[self.STORIcount] = {}
		self.runStats[self.STORIcount]["scores"] = []
		self.runStats[self.STORIcount]["scores"].append("-1")
		self.runStats[self.STORIcount]["runtime"] = float(-1)
	
	def __str__(self):                    #this method is called when print <Run> is used
		#print "Run data goes here!\n"
		print self.runName, ": "
		print self.__dict__
		return "\n"
			
	def StartRetrieval(self):
		subprocess.check_call(["perl", beginSTORIpath_t, self.runName, self.parentDir, self.taxaFile, str(self.windowSize), str(self.finalMaxFams)])
		#print "beginSTORI complete\n"
		self.SetRunIDs()
		#print "beginSTORI and SetRunIDs complete\n"
		return True
	
	
		'''
		ContinueRetrieval_t
		Once both iterators have finished, the STORI run controller finds the
	intersection of the results from each independent output, and provides the
	intersection set to each of a new pair of iterators. 
	'''
	def ContinueRetrieval_t(self):
		self.STORIcount += 1
		
		runNameA = self.runName + "a"
		runNameB = self.runName + "b"
		
		sourceFilesDirA = self.parentDir + "/" + runNameA
		sourceFilesDirB = self.parentDir + "/" + runNameB
		
		fileA = sourceFilesDirA + "/STORI_out_sc" + str(self.STORIcount-1) + "_" + runNameA + ".txt"
		fileB = sourceFilesDirB + "/STORI_out_sc" + str(self.STORIcount-1) + "_" + runNameB + ".txt"
	
		print "ContinueRetrieval_t cmd: " + " ".join(["perl", continueSTORIpath_t, fileA, fileB, runNameA, runNameB, self.parentDir, self.taxaFile, str(self.windowSize), str(self.finalMaxFams), str(self.STORIcount)])
		subprocess.check_call(["perl", continueSTORIpath_t, fileA, fileB, runNameA, runNameB, self.parentDir, self.taxaFile, str(self.windowSize), str(self.finalMaxFams), str(self.STORIcount)])
		
		self.runStats[self.STORIcount] = {}
		self.runStats[self.STORIcount]["scores"] = []
		self.runStats[self.STORIcount]["scores"].append("-1")
		self.runStats[self.STORIcount]["runtime"] = float(-1)
		
		self.SetRunIDs()
		return True
	
	def SetRunIDs(self):
		runIDfile = "".join([self.parentDir, "/tempRunIDfile.txt"])
		#print "attempting open ", runIDfile, "\n\n"
		
		m_runA = None
		m_runB = None
		beginTemp2 = open(runIDfile)
		for output in beginTemp2:
			if (m_runA is None):
				m_runA = (re.match('runA\s(\d+)\;', output))
			if (m_runB is None):
				m_runB = (re.match('runB\s(\d+)\;', output))
		beginTemp2.close()
		if (m_runA): self.runA_id = m_runA.group(1)
		if (m_runB): self.runB_id = m_runB.group(1)
		
		print "run ids ", self.runA_id, self.runB_id, "\n"
		
		self.runStats[self.STORIcount]["ids"] = (self.runA_id, self.runB_id)
	
	def CheckConvergence2(self, a, b, c):
		avg1 = int(100*(sum(map(float,self.runStats[a]["scores"]))/len(self.runStats[a]["scores"])))
		avg2 = int(100*(sum(map(float,self.runStats[b]["scores"]))/len(self.runStats[b]["scores"])))
		avg3 = int(100*(sum(map(float,self.runStats[c]["scores"]))/len(self.runStats[c]["scores"])))
		
		if ((avg1 > avg2) and (avg2 > avg3)):
			self.converged = False
		elif ((avg1 > 90) and (avg2 > 90) and (avg3 > 90)):
			if (abs(avg1-avg2)<4):
				self.converged = True
		
		#print jobStats "Check Convergence: $avg1 $avg2 $avg3 ";
	
		#if ($convergenceFlag==1) {
		#	print jobStats "convergence achieved.\n";
		#}
		#else {
		#	print jobStats "nope\n";
		#}
		

	
	
class RunSet():
	def __init__(self, arrayFile):
		self.arrayFile = arrayFile
		self.runHash = {}
		self.runArr = []
		self.ParseAndUpdateDone = True
		
		self.backupFileRWq_thread = BackupFileRWq(self)
		self.backupFileRWq_thread.start()
		self.backupFileRWq_thread.qLoadBackup()
		self.backupFileRWq_thread.qLoadParams()
		
	def __str__(self):                    #this method is called when print <RunSet> is used
		#return "RunSet data goes here!\n"
		print self.__dict__
		print "\n"
		for run in self.runHash:
			print self.runHash[run]
		return "\n"

	def LoadBackup(self):	#loads a specified CSV file into runArr
		self.runArr = []
		#print "runArr size: " + str(len(self.runArr)) + "\n"
		arrFile = None
		if (os.path.isfile(self.arrayFile)):
			arrFile = open(self.arrayFile)
			for line in arrFile:
				if (re.match('BEGIN\sPARAMETERS', line)): break
				line = line.rstrip()
				row = line.split(",")
				self.runArr.append(row)
			arrFile.close()
		
			#print "(end LoadBackup runArr)\n";
			#print "runArr size: " + str(len(self.runArr)) + "\n"
			#print self.runArr
			#print "\n"
		return True
		
	def LoadParams(self):
		self.runHash = {}
		arrFile = None
		readFlag = False
		tempParamsArr = []
		if (os.path.isfile(self.arrayFile)):
			arrFile = open(self.arrayFile)
			for line in arrFile:
				line = line.rstrip()
				if (re.match('BEGIN\sPARAMETERS', line)): 
					readFlag = True
					continue
				if readFlag:
					row = line.split(",")
					tempParamsArr.append(row)
			
			for tempRow in tempParamsArr:
				#print "tempRow in LoadParams: "
				#print tempRow
				#print "\n"
				runName = tempRow.pop(0)
				#print tempRow
				self.runHash[runName] = Run(runName, tempRow.pop(0), tempRow.pop(0), tempRow.pop(0), tempRow.pop(0), 0)


		
	def SaveBackup(self):			#output runHash and run_params to a CSV file at specified loc, overwriting any previous contents
		backupFile = open(self.arrayFile, "w")
		
		for runName in self.runHash.keys():
			STORIS_sorted = sorted(self.runHash[runName].runStats.keys())
			backupFile.write("".join([runName, ",", str(self.runHash[runName].runA_id), ",", str(self.runHash[runName].runB_id), ","]))
			for STORInum in STORIS_sorted:
				backupFile.write( \
					"".join( \
						[",".join(map(str, self.runHash[runName].runStats[STORInum]["scores"])), \
						",r", \
						str(self.runHash[runName].runStats[STORInum]["runtime"]), \
						","] \
					) \
				)
			if (self.runHash[runName].converged):
				backupFile.write("converged,")
			elif (self.runHash[runName].paused):
				backupFile.write("paused,")
			backupFile.write("\n")
		
		backupFile.write("BEGIN PARAMETERS\n")
		
		for runName in self.runHash.keys():
			backupFile.write("".join([runName, ",", self.runHash[runName].parentDir, ",", self.runHash[runName].taxaFile, ",", str(self.runHash[runName].windowSize), ",", str(self.runHash[runName].finalMaxFams), "\n"]))
		
		backupFile.close()
		return True
		
	
	def StartRetrieval_t(self, runName, parentDir, taxaFile, windowSize, finalMaxFams, STORIcount):
		print "starting runName=", runName, " parentDir=", parentDir, " taxaFile=", taxaFile, " windowSize=", windowSize, " finalMaxFams=", finalMaxFams, " STORIcount=", STORIcount, "\n"
		if not (runName in self.runHash):
			self.runHash[runName] = Run(runName, parentDir, taxaFile, windowSize, finalMaxFams, STORIcount)
			methodDone = False
			methodDone = self.backupFileRWq_thread.qStartRetrieval(runName)   #we don't want a StartRetrieval to collide with a ContinueRetrieval since they both write the temp pid file
			#print "waiting on qStartRetrieval\n"
			while (not methodDone):
				time.sleep(1)
		else:
			print runName + " already exists!\n"
		return True
		
	def StopRun(self, runName):
		if runName in self.runHash:
			print "stopping ", runName, " and deleting it from runHash. Runfiles are still present.\n"
			try: subprocess.check_call(["kill", str(self.runHash[runName].runA_id)])
			except subprocess.CalledProcessError: print "PID " + str(self.runHash[runName].runA_id) + " no longer exists\n"
			#print "first kill attempted\n"
			try: subprocess.check_call(["kill", str(self.runHash[runName].runB_id)])
			except subprocess.CalledProcessError: print "PID " + str(self.runHash[runName].runB_id) + " no longer exists\n"
			#print "second kill attempted\n"
			del self.runHash[runName]
			#print "Run deleted from runHash\n"
		else:
			print runName + " not found.\n"
		return True	

		
	def PauseRun(self, runName):
		if runName in self.runHash:
			print "pausing ", runName, " killing PIDs, but it remains in runHash. Runfiles are still present.\n"
			try: subprocess.check_call(["kill", str(self.runHash[runName].runA_id)])
			except subprocess.CalledProcessError: print "PID " + str(self.runHash[runName].runA_id) + " no longer exists\n"
			
			try: subprocess.check_call(["kill", str(self.runHash[runName].runB_id)])
			except subprocess.CalledProcessError: print "PID " + str(self.runHash[runName].runB_id) + " no longer exists\n"
			
			self.runHash[runName].paused = True
			#print "scores for " + runName + ": "
			#print self.runHash[runName].runStats[self.runHash[runName].STORIcount]["scores"]
			#print "\n"
		else:
			print runName + " not found.\n"
		return True
					
	def ParseAndUpdate_caller(self):
		while True:
			self.ParseAndUpdateDone = False
			#print "beginning self.ParseAndUpdate() qSaveBackup() qLoadBackup() self.ParseAndUpdateDone=" + str(self.ParseAndUpdateDone) + "\n"
			self.ParseAndUpdateDone = self.ParseAndUpdate()
			
			self.backupFileRWq_thread.qSaveBackup()
			self.backupFileRWq_thread.qLoadBackup() #need a load here b/c need to have correct pids in runArr
			#print "finished self.ParseAndUpdate() qSaveBackup() qLoadBackup() self.ParseAndUpdateDone=" + str(self.ParseAndUpdateDone) + "\n"
			
			time.sleep(15)
	
	def ParseAndUpdate(self):
		#print "parsing and updating yo\n"
		for row in self.runArr:
			tempRow = list(row)
			tempRow.pop() # remove blank last element
			runName = tempRow.pop(0)
			runid_A = tempRow.pop(0)
			runid_B = tempRow.pop(0)
			score = ""
			STORIcount = int(0)
			
			while tempRow:
				scores = []
				#print "[355] score=" + str(score) + " tempRow=" + str(tempRow) + "\n"
				while (not(re.match('r', score))):
					if (re.match('-*\d+', score)):
						scores.append(score)
					score = tempRow.pop(0)
					#print "[360] score=" + str(score) + " tempRow=" + str(tempRow) + "\n"
				
				#print "ParseAndUpdate here are scores: "
				#print scores
				
				self.runHash[runName].runStats[STORIcount] = {}
				self.runHash[runName].runStats[STORIcount]["scores"] = list(scores)
				self.runHash[runName].STORIcount = STORIcount
				
				score = score.rstrip()
				rtime = re.match('r(.+)', score).group(1)
				
				#print "<ParseAndUpdate> self.runHash[" + runName + "].runStats[" + str(STORIcount) + "][\"runtime\"] = " + rtime + "\n"
				self.runHash[runName].runStats[STORIcount]["runtime"] = float(rtime)
				#if (float(rtime) > 0): score = tempRow.pop(0)		#to exclude the case where the runtime is -1, and therefore tempRow is empty 
				
				if tempRow:		#if not empty
					score = tempRow.pop(0)
					if (re.match('converged', score)):
						self.runHash[runName].converged = True
					if (re.match('paused', score)):
						self.runHash[runName].paused = True
				
				STORIcount+=1
			
			STORIcount-=1
			
			self.runHash[runName].runA_id = int(runid_A)
			self.runHash[runName].runB_id = int(runid_B)
			
			'''
			Periodically, the STORI controller compares the results of the iterators, and
			assigns the pair a score reflecting how similar their orthology predictions are
			to one another. 
			'''
			
		for runName in self.runHash.keys():			#cycle through the not-converged runs in runHash, updating scores and runtimes, and start a new run if prev are finished and not converged
			if ((self.runHash[runName].converged == False) and (self.runHash[runName].paused == False)):
				#print "checking non-converged non-paused runs\n"
				sourceFilesDirA = self.runHash[runName].parentDir + "/" + runName + "a"
				sourceFilesDirB = self.runHash[runName].parentDir + "/" + runName + "b"
				
				STORIcount = self.runHash[runName].STORIcount

				fileA = sourceFilesDirA + "/STORI_out_sc" + str(STORIcount) + "_" + runName + "a.txt"
				fileB = sourceFilesDirB + "/STORI_out_sc" + str(STORIcount) + "_" + runName + "b.txt"
				
				output = subprocess.check_output(["perl", checkSTORIpath, fileA, fileB, "-q"])
				#print str(output) + "\n"
				runtime_a = re.search('runtime_a:\s(.+)\n', output).group(1)
				runtime_b = re.search('runtime_b:\s(.+)\n', output).group(1)
				
				runtime_a_f = float(-1)
				runtime_b_f = float(-1)
				
				if not(re.match('^-1$', runtime_a)):
					rtime_m = re.match('(\d+):(\d+):(\d+)', runtime_a)
					hours = float(rtime_m.group(1))
					hours_mins = (float(rtime_m.group(2))/60)
					hours_secs = (float(rtime_m.group(3))/3600)
					runtime_a_f = hours + hours_mins + hours_secs
				
				if not(re.match('^-1$', runtime_b)):
					rtime_m = re.match('(\d+):(\d+):(\d+)', runtime_b)
					hours = float(rtime_m.group(1))
					hours_mins = (float(rtime_m.group(2))/60)
					hours_secs = (float(rtime_m.group(3))/3600)
					runtime_b_f = hours + hours_mins + hours_secs
	
				if ((runtime_a_f > 0) and (runtime_b_f > 0)):
					cpuTime = runtime_a_f + runtime_b_f	
					self.runHash[runName].runStats[STORIcount]["runtime"] = cpuTime
					#print "runHash{$runName}{$STORIcount}{runtime}= $cpuTime \n";
	
				#runtimes have been established. now time to update scores
				#print "(end of ParseAndUpdate) output = $output \n";
				#print "score match=" + re.search('score\s(.+)\s\n',output).group(1) + "END\n"
				scores = re.search('score\s(.+)\s\n',output).group(1).split(" ")
				#print "scores= " + str(scores) + "\n"
				self.runHash[runName].runStats[STORIcount]["scores"] = map(float, scores)
				
				#print "runHash{$runName}{$STORIcount}{scores} = [" . join(" ",@scores) . "]\n";
				
				STORIfinished = re.search('bothfinished', output)
				
				if (((runtime_a_f > maxRuntime) and (runtime_b_f > maxRuntime)) or (STORIfinished)):
					print "either maxRuntime exceeded or STORI.pl bothfinished\n"
					print "STORIfinished=" + str(bool(STORIfinished)) + "\n"
					print "maxRuntime=" + str(maxRuntime) + "  runtime_a_f=" + str(runtime_a_f) + "  runtime_b_f=" + str(runtime_b_f) + "\n"
					try: subprocess.check_call(["kill", str(self.runHash[runName].runA_id)])
					except subprocess.CalledProcessError: print "PID " + str(self.runHash[runName].runA_id) + " no longer exists\n"
					
					try: subprocess.check_call(["kill", str(self.runHash[runName].runB_id)])
					except subprocess.CalledProcessError: print "PID " + str(self.runHash[runName].runB_id) + " no longer exists\n"
					
					methodDone = False
					if (STORIcount < 2):
						print "waiting on qContinueRetrieval\n"
						methodDone = self.backupFileRWq_thread.qContinueRetrieval(runName)   #we don't want a StartRetrieval to collide with a ContinueRetrieval since they both write the temp pid file
						while (not methodDone):
							time.sleep(1)
					elif (STORIcount >= 2):
						a = STORIcount
						b = (a - 1)
						c = (b - 1)
						print "checking convergence...\n"
						'''
						Eventually the similarity score stabilizes in the range of
						90-100%. When the controller detects this stabilization, it labels the run as
						converged and stops further iteration.
						'''
						self.runHash[runName].CheckConvergence2(a,b,c)
						if not(self.runHash[runName].converged):
							#tempScore = (sum(self.runHash[runName].runStats[STORIcount]["scores"])/len(self.runHash[runName].runStats[STORIcount]["scores"]))
							methodDone = self.backupFileRWq_thread.qContinueRetrieval(runName)
							#print "waiting on qContinueRetrieval\n"
							while (not methodDone):
								time.sleep(1)
						elif self.runHash[runName].converged:
							print "removing temp blast db from source dir\n"
							delPathA = self.runHash[runName].parentDir + "/" + runName + "a/tempBlastDB"
							print "cmd is: rm -r " + delPathA + "\n"
							subprocess.check_call(["rm","-r",delPathA])
							delPathB = self.runHash[runName].parentDir + "/" + runName + "b/tempBlastDB"
							print "cmd is: rm -r " + delPathB + "\n"
							subprocess.check_call(["rm","-r",delPathB])
							#remove the database to conserve disk space
						
		return True

class BackupFileRWq(threading.Thread):		#inherits threading.Thread class
											#not just for backup file - also for pid storage file access
		def __init__(self, RunSetObj):
			self.q = Queue.Queue()
			self.RunSetObj = RunSetObj
			threading.Thread.__init__(self)
			self.daemon = True
			self.methodDone = False
			
		def run(self):
			while True:
				while (not self.RunSetObj.ParseAndUpdateDone):
					time.sleep(1)
				time.sleep(1)
				#print "in BackupFileRWq loop\n"
				val = self.q.get()
				#print "val=" + str(val) + "\n"
				if (val == 1): 
					#print "q1 LoadParams()\n"
					self.RunSetObj.LoadParams()
				elif (val == 2): 
					#print "q2 LoadBackup()\n"
					self.methodDone = self.RunSetObj.LoadBackup()
				elif (val == 3): 
					#print "q3 SaveBackup()\n"
					self.RunSetObj.SaveBackup()
				elif (val == 4): 
					#print "q4 StartRetrieval()\n"
					runName = self.q.get()
					self.methodDone = self.RunSetObj.runHash[runName].StartRetrieval()
					self.RunSetObj.SaveBackup()
					self.RunSetObj.LoadBackup()
				elif (val == 5):
					#print "q5 ContinueRetrieval_t()\n"
					runName = self.q.get()
					self.methodDone = self.RunSetObj.runHash[runName].ContinueRetrieval_t()
					self.RunSetObj.SaveBackup()
					self.RunSetObj.LoadBackup()
				elif (val == 6):
					#print "q6 StopRun()\n"
					runName = self.q.get()
					self.methodDone = self.RunSetObj.StopRun(runName)
					#print "q6 StopRun() methodDone= " + str(self.methodDone) + "\n"
					self.RunSetObj.SaveBackup()
					#print "q6 StopRun()Saved backup\n"
					self.RunSetObj.LoadBackup()
					#print "q6 StopRun loaded backup()\n"
				elif (val == 7):
					#print "q7 PauseRun()\n"
					runName = self.q.get()
					self.methodDone = self.RunSetObj.PauseRun(runName)
					self.RunSetObj.SaveBackup()
					self.RunSetObj.LoadBackup()

		def qLoadParams(self):
			self.q.put(1)
		def qLoadBackup(self):
			self.q.put(2)
			while (not self.methodDone):
				time.sleep(1)
			self.methodDone = False
			return True
		def qSaveBackup(self):
			self.q.put(3)
		def qStartRetrieval(self, runName):
			self.q.put(4)
			self.q.put(runName)
			while (not self.methodDone):
				time.sleep(2)
			self.methodDone = False
			return True
				
		def qContinueRetrieval(self, runName):
			self.q.put(5)
			self.q.put(runName)
			#print "entered qContinueRetrieval after put, self.methodDone=" + str(self.methodDone) + " runName=" + runName + " qsize=" + str(self.q.qsize()) + "\n"
			
			while (not self.methodDone):
				time.sleep(5)
				#print "qsize=" + str(self.q.qsize()) + " qContinueRetrieval sleeping 5 sec\n"
			self.methodDone = False
			#print "qContinueRetrieval finished sleeping, self.methodDone=" + str(self.methodDone) + "\n"
			return True
			
		def qStopRun(self, runName):
			self.q.put(6)
			self.q.put(runName)
			while (not self.methodDone):
				time.sleep(2)
			self.methodDone = False
			return True
		def qPauseRun(self, runName):
			self.q.put(7)
			self.q.put(runName)
			while (not self.methodDone):
				time.sleep(2)
			self.methodDone = False
			return True



runs1 = RunSet(arrayFile)

print  """\n\t\tWelcome to STORI Control!
		  \n\t\tMulti-node STORI (c)2012 JG Stern & EA Gaucher
\t\tSingle-node STORI (c)2015 JG Stern
	   \n\t\t\tCommands:
	   \t\t\tstart <run-name> <scratch/dir> <taxa file> <windowSize> <finalMaxFams>
	   \t\t\tstop <run-name>
	   \t\t\tpause <run-name>
	   \t\t\texit\n"""
	   
#!!!!!!!
#starts the runs1.ParseAndUpdate_caller() thread. will run in background. every 5 seconds or so
ParseAndUpdateThread = threading.Thread(target=runs1.ParseAndUpdate_caller)
ParseAndUpdateThread.daemon = True
ParseAndUpdateThread.start()


while True:
	time.sleep(.1)
	#print "in the main loop\n"
	methodDone = False
	inputStr = raw_input("STORI>")
	
	m_show = (re.match('show\s(.+)', inputStr))
	m_start = (re.match('start\s(?P<runName>.+)\s(?P<parentDir>.+)\s(?P<taxaFile>.+)\s(?P<windowSize>.+)\s(?P<finalMaxFams>.+)', inputStr))
	m_stop = (re.match('stop\s(?P<runName>.+)', inputStr))
	m_pause = (re.match('pause\s(?P<runName>.+)', inputStr))
	m_exit = (re.match('exit', inputStr))
	m_newline = (re.match('^$', inputStr))
	
	if m_show:
		if (re.match('groups', m_show.group(1))):
			print "showing the groups\n"
		else:
			print "showing the runs\n"
			print runs1.runHash.__dict__

	elif m_start:
		methodDone = runs1.StartRetrieval_t(m_start.group('runName'), m_start.group('parentDir'), m_start.group('taxaFile'), m_start.group('windowSize'), m_start.group('finalMaxFams'), 0)
		while (not methodDone):
			time.sleep(1)
			#print "StartRetrieval_t methodDone= " + str(methodDone) + "\n"
		#print "StartRetrieval_t methodDone= " + str(methodDone) + "\n"
	elif m_stop:
		runs1.backupFileRWq_thread.qStopRun(m_stop.group('runName'))
	elif m_pause:
		runs1.backupFileRWq_thread.qPauseRun(m_pause.group('runName'))
	elif m_exit:
		print "bye!\n"
		break
	elif m_newline:
		print ""
	else:
		print "Sorry, I don't understand.\n"
		print  """
		   \n\t\t\tCommands:
		   \t\t\tstart <run-name> <scratch/dir> <taxa file> <windowSize> <finalMaxFams>
		   \t\t\tstop <run-name>
		   \t\t\tpause <run-name>
		   \t\t\texit\n"""
	
	#print runs1

#print "exited main input loop\n"
