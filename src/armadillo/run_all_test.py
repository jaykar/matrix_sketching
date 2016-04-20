import glob
import os
import subprocess

list_files = glob.glob('test/*')
for file_name in list_files:
    command = "make" 
    print command, file_name
    (sd, se) = subprocess.Popen([command], stdout=subprocess.PIPE, cwd=file_name).communicate()
    print "stdout, stderr:"
    print sd, se


    #hacky path to valgrind, regular valgrind is broken on my computer and I didn't want to change that path
    #so I pointed to a dev version
    command2 = '/usr/local/Cellar/valgrind/HEAD/bin/valgrind --leak-check=yes ./' + file_name + '/sample1_unittest'
    print command2
    (sd, se) = subprocess.Popen([command2], stdout=subprocess.PIPE, shell=True).communicate()
    print "stdout, stderr"
    print sd, se
