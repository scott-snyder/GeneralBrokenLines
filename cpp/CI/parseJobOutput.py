import sys
from argparse import ArgumentParser 
import testDefs 

def getLines(theFileName: str):
    with open(theFileName, "r") as fin:
        # we do not validate lines showing run times, as those are not expected to be 
        # 1:1 reproducible. 
        li = [l.rstrip("\n") for l in fin.readlines() if not "Time elapsed" in l]
        return li 

def compare(testName: str, testDir: str, refDir: str):
    print (f"Validating test {testName}")
    refContent = getLines(f"{refDir}/{testDefs.logFileName(testName, True)}")
    testContent = getLines(f"{testDir}/{testDefs.logFileName(testName, False)}")
    
    # for now, we perform a very strict check for full identity. 
    # Might want to relax in the future.

    stat = 0  
    if len(refContent) != len(testContent):
        print (f"Files differ in line count: Ref has {len(refContent)}, test has {len(testContent)}")
        stat = 1
    for i, refToken in enumerate (refContent):
        if refToken not  in testContent: 
            print (f"Line {i} of the ref file not found in test:\n   {refToken}")
            stat = 1 
    for i, testToken in enumerate (testContent):
        if testToken not  in refContent: 
            print (f"Line {i} of the test file not found in ref:\n   {testToken}")
            stat = 1 
        
    if (stat == 0):
        print ("    Ref and test in exact agreement")
        return True 
    else:
        print ("    Ref and test NOT in exact agreement.")
        return False 

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("-t", dest="tests",default=testDefs.listOfTests,nargs="+", help="specify tests to compare")
    parser.add_argument("testDir",type=str, help="Location of the test files")
    parser.add_argument("refDir",type=str, help="Location of the ref files")

    args = parser.parse_args()

    stat = 0
    nOK = 0
    nTests = len(args.tests)
    for test in args.tests:
        if compare(test, args.testDir, args.refDir):
            nOK+=1
        else:
            stat = 1 
    print (f"=== {nOK} out of {nTests} comparisons succeeded ===")
    sys.exit(stat)
