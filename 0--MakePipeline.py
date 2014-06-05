#!/usr/bin/env python

import getopt,sys
import scriptMaker

def usage():
    print "Outputs all pipeline scripts given a configuration file."
    print "Usage: 0--MakePipeline.py [-h --help] [-c --conf <configuration file>]"

try:
    opts, args = getopt.getopt(sys.argv[1:], "hc:", ["help", "conf"])
except getopt.GetoptError:
    usage()
    sys.exit(2)

## Default settings
conf=None
####################

for o, a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit()
    if o in ("-c", "--conf"):
        conf=a

#####################        
if conf:
    print "Starting..."

    try:
        f=open(conf,"rU")
    except IOError, e:
        print "File not found: [", conf, "]"
        sys.exit(2)

    param={}
    for i in f.readlines():
        tmp=i.strip().split("\t")
        param[tmp[0]]=tmp[1]
    f.close()

    scriptMaker.FastQC(param)
    scriptMaker.MappingAndPreProcessing(param)
    scriptMaker.QualityControl(param)
    scriptMaker.HaplotypeCaller(param)

    print "Done."

####################    
else:
    usage()
    sys.exit()

