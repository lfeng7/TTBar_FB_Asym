import glob
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-d','--directory', metavar='F', type='string', action='store',dest='directory',help='') ## Sets which files to run on
(options, args) = parser.parse_args()

FILES = options.directory + '/*.root'
print FILES
files = glob.glob(FILES)
for fname in files :
    cmd = 'echo '+fname+' >> input_files_list.txt'
    print cmd
    os.system(cmd)
