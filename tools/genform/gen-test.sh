genform ms=test-data/CreatineMs.txt out=test-data/CreatineMs.out exist out=test-data/CreatineMs.out acc=20 rej=30 

genform ms=test-data/CreatineMs.txt msms=test-data/CreatineMsMs.txt exist oms=test-data/CreatineMsMs-oms.out omsms=test-data/CreatineMsMs-omsms.out oclean=test-data/CreatineMsMs-oclean.out out=test-data/CreatineMsMs.out acc=20 rej=30 analyze loss intens dbe cm pc sc
