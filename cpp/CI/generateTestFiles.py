import subprocess 
import sys
from argparse import ArgumentParser

import testDefs 

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("-r",action="store_true",dest="refFiles",default=False, help="produce ref instead of test files")
    parser.add_argument("-t", dest="tests",default=testDefs.listOfTests,nargs="+", help="specify tests to run")
    parser.add_argument("-n", dest="nTracks",default=100, help="Number of tracks to fit per test")
    args = parser.parse_args()
    
    print (f" Will run the following tests: {args.tests}")
    
    if (args.refFiles): 
        print ("Will generate ref files")
    status = 0 
    for test in args.tests:
        stat, log = subprocess.getstatusoutput(f"{test} {args.nTracks}") 
        if stat != 0:
            print (f" FAILED to run test {test} {args.nTracks}")
            status = 1 
        else: 
            print (f"Test {test} successfully run")
            with open(testDefs.logFileName(test, args.refFiles) ,"w") as fout: 
                fout.write(log)
    sys.exit( status) 
        

    
