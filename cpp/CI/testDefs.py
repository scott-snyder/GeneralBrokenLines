listOfTests = ["GBL_example1",
                       "GBL_example2",
                       "GBL_example3",
                       "GBL_example4", 
                       "GBL_exampleComposedGeo", 
                       "GBL_exampleComposedKin", 
                       "GBL_exampleDc", 
                       "GBL_exampleSit"
                       ]

def logFileName(testName: str, isRef: bool = False)-> str:
    filePrefix="test."
    if (isRef): 
        filePrefix="ref."

    return f"{filePrefix}{testName}.log"