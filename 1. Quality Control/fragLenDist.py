import h5py
import numpy as np
import pandas as pd
#from collections import defaultdict

snap_fp = "/home/chenyin/project/mouseOlfactory/mmob.snap"
snap = h5py.File(snap_fp, "r")
fragLen = np.array(snap["FM/fragLen"])
fragLen = pd.Series(fragLen)
counts = fragLen.value_counts()
counts = dict(counts)

opstr = '\n'.join(['%s\t%s'%(k,v) for k,v in counts.items()])

with open("fragLenStat.txt", "w") as writer:
    writer.write(opstr)
    
    
