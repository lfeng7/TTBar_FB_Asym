import os
import glob

#first get the list of all the output files
outputfilelist = glob.glob('output_JID_*.log')
#make a list of the failed job numbers
failedjobnumbers = []
#start by checking each output file to see if the job failed
for outputfile in outputfilelist :
	jobend = os.popen('tail -n 5 '+outputfile+' | head -n 1').read()
	if not jobend.startswith('Count at ') :
		if jobend.find('PROBLEM IN FIT: ierflag = 4')!=-1 :
			continue
		if jobend.find(' Xrd: XrdClientMessage::ReadRaw: Failed to read header')!=-1 :
			continue
		if jobend.find('Xrd: CheckErrorStatus: Server [cmseos.fnal.gov:')!=-1 :
			continue
		jobnumber = outputfile.rstrip('.log').split('_')[len(outputfile.rstrip('.log').split('_'))-1]
		print 'Job '+jobnumber+' failed with last line "'+jobend.rstrip('\n')+'"'
		failedjobnumbers.append(int(jobnumber))
#now check the file sizes to find any that are abnormally small
rootfilelist = glob.glob('*_tree.root')
totalsize = 0.
for rootfile in rootfilelist :
	totalsize+=os.path.getsize(rootfile)
expected_contribution = totalsize/len(rootfilelist)
for rootfile in rootfilelist :
	filesize = os.path.getsize(rootfile)
	if filesize/expected_contribution<0.95 :
		print 'File '+rootfile+' is too small, its size is '+str(filesize)+' bytes, contributing '+str(filesize/expected_contribution)+' of its expectation'
		jobnumber = int(rootfile.rstrip('_tree.root').split('_')[len(rootfile.rstrip('_tree.root').split('_'))-1])
		if not 'Run2012' in rootfile :
			jobnumber *= 5
			if rootfile.find('JES_up')!=-1 :
				jobnumber+=1
			elif rootfile.find('JES_down')!=-1 :
				jobnumber+=2
			elif rootfile.find('JER_up')!=-1 :
				jobnumber+=3
			elif rootfile.find('JER_down')!=-1 :
				jobnumber+=4
		if not jobnumber in failedjobnumbers :
			failedjobnumbers.append(jobnumber)
#now check if any jobs have output but no returned file
if len(outputfilelist) > len(rootfilelist) :
	rootfilejobnumbers = []
	for rootfile in rootfilelist :
		rootfilejobnumber = int(rootfile.rstrip('_tree.root').split('_')[len(rootfile.rstrip('_tree.root').split('_'))-1])
		if not 'Run2012' in rootfile :
			rootfilejobnumber *= 5
			if rootfile.find('JES_up')!=-1 :
				rootfilejobnumber+=1
			elif rootfile.find('JES_down')!=-1 :
				rootfilejobnumber+=2
			elif rootfile.find('JER_up')!=-1 :
				rootfilejobnumber+=3
			elif rootfile.find('JER_down')!=-1 :
				rootfilejobnumber+=4
		rootfilejobnumbers.append(rootfilejobnumber)
	outputfilejobnumbers = []
	for outputfile in outputfilelist :
		outputfilejobnumber = int(outputfile.rstrip('.log').split('_')[len(outputfile.rstrip('.log').split('_'))-1])
		outputfilejobnumbers.append(outputfilejobnumber)
	for n in outputfilejobnumbers :
		if n not in rootfilejobnumbers :
			print 'job number '+str(n)+' has no outputted root file.'
			failedjobnumbers.append(n)
#sort the list of failed job numbers
failedjobnumbers.sort()
#open the list of all the jobs and add the failed ones to the new file
linecount = 0
if not os.path.isfile('ana.listOfJobs_all') :
	print 'TOTAL LIST OF JOBS DOES NOT EXIST YET, COPYING CURRENT LIST OF JOBS!!'
	os.system('mv ana.listOfJobs ana.listOfJobs_all')
os.system('rm -rf ana.listOfJobs')
joblist = open('ana.listOfJobs_all','r')
for job in joblist.readlines() :
	if linecount in failedjobnumbers :
		os.system('echo "'+job+'" >> ana.listOfJobs')
	linecount+=1
print 'Total new list of jobs: '
os.system('cat ana.listOfJobs')
os.system('bash cleanup.bash')
