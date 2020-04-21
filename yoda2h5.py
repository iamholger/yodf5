import yoda, sys
import h5py
import numpy as np


def createDatasets(f, binids, variations , compression=4):
    nbins=len(binids)
    nvars=len(variations)

    floats = [
            "sumw",
            "sumw2",
            "sumwx",
            "sumwx2",
            "sumwy",
            "sumwy2",
            "sumwxy",
            "numEntries",
            "xval",
            "xerr-",
            "xerr+",
            "yval",
            "yerr-",
            "yerr+",
            "xmin",
            "xmax",
            "ymin",
            "ymax"
            ]

    for df in floats: f.create_dataset(df, (nbins,nvars), maxshape=(None,None), dtype='f' , chunks=True, compression=compression)

    # Lookups
    f.create_group("Histo1D")
    f.create_group("Histo2D")
    f.create_group("Profile1D")
    f.create_group("Counter")
    f.create_group("Scatter1D")
    f.create_group("Scatter2D")

    # This is the one that works well with hdf5 when reading std::string
    dt = h5py.special_dtype(vlen=str)
    import numpy as np
    f.create_dataset("binids",     data=np.array(binids,     dtype=dt))
    f.create_dataset("variations", data=np.array(variations, dtype=dt))


def dbn0ToArray(dbn):
    return np.array([dbn.sumW(), dbn.sumW2(), dbn.numEntries()])

def dbn1ToArray(dbn):
    """
    The try except block deals with the underflow things not hving xmin, xmax
    """
    try:
        return np.array([dbn.sumW(), dbn.sumW2(), dbn.sumWX(), dbn.sumWX2(), dbn.numEntries(), dbn.xMin(), dbn.xMax()])
    except:
        return np.array([dbn.sumW(), dbn.sumW2(), dbn.sumWX(), dbn.sumWX2(), dbn.numEntries(), 0, 0])

def H2dbn2ToArray(dbn):
    try:
        return np.array([dbn.sumW(), dbn.sumW2(), dbn.sumWX(), dbn.sumWX2(), dbn.sumWY(), dbn.sumWY2(), dbn.sumWXY(), dbn.numEntries(), dbn.xMin(), dbn.xMax(), dbn.yMin(), dbn.yMax()])
    except:
        return np.array([dbn.sumW(), dbn.sumW2(), dbn.sumWX(), dbn.sumWX2(), dbn.sumWY(), dbn.sumWY2(), dbn.sumWXY(), dbn.numEntries(), 0, 0, 0, 0])

def dbn2ToArray(dbn):
    try:
        return np.array([dbn.sumW(), dbn.sumW2(), dbn.sumWX(), dbn.sumWX2(), dbn.sumWY(), dbn.sumWY2(), dbn.numEntries(), dbn.xMin(), dbn.xMax()])
    except:
        return np.array([dbn.sumW(), dbn.sumW2(), dbn.sumWX(), dbn.sumWX2(), dbn.sumWY(), dbn.sumWY2(), dbn.numEntries(), 0, 0])

def point2DToArray(pnt):
    return np.array([pnt.val(1), pnt.errMinus(1), pnt.errPlus(1), pnt.val(2), pnt.errMinus(2), pnt.errPlus(2)])

def point1DToArray(pnt):
    return np.array([pnt.val(1), pnt.errMinus(1), pnt.errPlus(1)])

def fillDatasets(f, binids, variations, ddict, hname):
    _hname=hname.replace("/","|")

    binIdx = [num for num, binid in enumerate(binids) if binid.startswith("{}#".format(_hname))]
    if len(binIdx) ==0:
        print("Warning, no matching binid for {} --- is this one of the raw ratios maybe???".format(hname))
        return

    hids = [hname]
    for v in variations[1:]:
        hids.append("{}[{}]".format(hname, v))

    if ddict[hname].type()=='Histo1D':
        nFields=7
        fdbn = dbn1ToArray
        f.create_dataset("Histo1D/{}".format(_hname), data=binIdx , chunks=True)
    elif ddict[hname].type()=='Histo2D':
        nFields=12
        fdbn = H2dbn2ToArray
        f.create_dataset("Histo2D/{}".format(_hname), data=binIdx , chunks=True)
    elif ddict[hname].type()=='Profile1D':
        fdbn = dbn2ToArray
        nFields=9
        f.create_dataset("Profile1D/{}".format(_hname), data=binIdx , chunks=True)
    elif ddict[hname].type()=='Scatter2D':
        fdbn = point2DToArray
        nFields=6
        f.create_dataset("Scatter2D/{}".format(_hname), data=binIdx , chunks=True)
    elif ddict[hname].type()=='Scatter1D':
        fdbn = point1DToArray
        nFields=3
        f.create_dataset("Scatter1D/{}".format(_hname), data=binIdx , chunks=True)
    elif ddict[hname].type()=='Counter':
        nFields=3
        f.create_dataset("Counter/{}".format(_hname), data=binIdx , chunks=True)
    else:
        raise Exception("type {} Not implemented".format(ddict[hname].type()))


    # This covers all non scatters
    temp = np.zeros((len(binIdx), len(variations), nFields))
    #sumw	 sumw2	 sumwx	 sumwx2	 sumwy	 sumwy2	 numEntries

    # Iterate over variations
    for col, hn in enumerate(hids):
        # Iterate over bins
        H=ddict[hn]

        if H.type() == "Counter":
            temp[0][col] = np.array([H.sumW(), H.sumW2(), H.numEntries()])

        # Things with under/overflow first
        elif H.type() not in ["Scatter1D", "Scatter2D", "Histo2D"]:
            temp[0][col] = fdbn(H.totalDbn())
            temp[1][col] = fdbn(H.overflow())
            temp[2][col] = fdbn(H.underflow())
            for i in range(len(binIdx)-3):
                temp[3+i][col] = fdbn(H.bin(i))
        elif H.type() =="Histo2D":
            temp[0][col] = fdbn(H.totalDbn())
            temp[1][col] = 0.0 # Future proofing
            temp[2][col] = 0.0 # 
            for i in range(len(binIdx)-3):
                temp[3+i][col] = fdbn(H.bin(i))
        else:
            for i in range(len(binIdx)):
                temp[i][col] = fdbn(H.point(i))

    if ddict[hname].type()=='Histo1D':
        f["sumw"][      binIdx[0]:binIdx[-1]+1,:] = temp[:,:,0]
        f["sumw2"][     binIdx[0]:binIdx[-1]+1,:] = temp[:,:,1]
        f["sumwx"][     binIdx[0]:binIdx[-1]+1,:] = temp[:,:,2]
        f["sumwx2"][    binIdx[0]:binIdx[-1]+1,:] = temp[:,:,3]
        f["numEntries"][binIdx[0]:binIdx[-1]+1,:] = temp[:,:,4]
        f["xmin"][      binIdx[0]:binIdx[-1]+1,:] = temp[:,:,5]
        f["xmax"][      binIdx[0]:binIdx[-1]+1,:] = temp[:,:,6]

    elif ddict[hname].type()=='Histo2D':
        f["sumw"][      binIdx[0]:binIdx[-1]+1,:] = temp[:,:,0]
        f["sumw2"][     binIdx[0]:binIdx[-1]+1,:] = temp[:,:,1]
        f["sumwx"][     binIdx[0]:binIdx[-1]+1,:] = temp[:,:,2]
        f["sumwx2"][    binIdx[0]:binIdx[-1]+1,:] = temp[:,:,3]
        f["sumwy"][     binIdx[0]:binIdx[-1]+1,:] = temp[:,:,4]
        f["sumwy2"][    binIdx[0]:binIdx[-1]+1,:] = temp[:,:,5]
        f["sumwxy"][    binIdx[0]:binIdx[-1]+1,:] = temp[:,:,6]
        f["numEntries"][binIdx[0]:binIdx[-1]+1,:] = temp[:,:,7]
        f["xmin"][      binIdx[0]:binIdx[-1]+1,:] = temp[:,:,8]
        f["xmax"][      binIdx[0]:binIdx[-1]+1,:] = temp[:,:,9]
        f["ymin"][      binIdx[0]:binIdx[-1]+1,:] = temp[:,:,10]
        f["ymax"][      binIdx[0]:binIdx[-1]+1,:] = temp[:,:,11]

    elif ddict[hname].type()=='Profile1D':
        f["sumw"][      binIdx[0]:binIdx[-1]+1,:] = temp[:,:,0]
        f["sumw2"][     binIdx[0]:binIdx[-1]+1,:] = temp[:,:,1]
        f["sumwx"][     binIdx[0]:binIdx[-1]+1,:] = temp[:,:,2]
        f["sumwx2"][    binIdx[0]:binIdx[-1]+1,:] = temp[:,:,3]
        f["sumwy"][     binIdx[0]:binIdx[-1]+1,:] = temp[:,:,4]
        f["sumwy2"][    binIdx[0]:binIdx[-1]+1,:] = temp[:,:,5]
        f["numEntries"][binIdx[0]:binIdx[-1]+1,:] = temp[:,:,6]
        f["xmin"][      binIdx[0]:binIdx[-1]+1,:] = temp[:,:,7]
        f["xmax"][      binIdx[0]:binIdx[-1]+1,:] = temp[:,:,8]

    elif ddict[hname].type()=='Scatter1D':
        f["xval"][ binIdx[0]:binIdx[-1]+1,:] = temp[:,:,0]
        f["xerr-"][binIdx[0]:binIdx[-1]+1,:] = temp[:,:,1]
        f["xerr+"][binIdx[0]:binIdx[-1]+1,:] = temp[:,:,2]

    elif ddict[hname].type()=='Scatter2D':
        f["xval"][ binIdx[0]:binIdx[-1]+1,:] = temp[:,:,0]
        f["xerr-"][binIdx[0]:binIdx[-1]+1,:] = temp[:,:,1]
        f["xerr+"][binIdx[0]:binIdx[-1]+1,:] = temp[:,:,2]
        f["yval"][ binIdx[0]:binIdx[-1]+1,:] = temp[:,:,3]
        f["yerr-"][binIdx[0]:binIdx[-1]+1,:] = temp[:,:,4]
        f["yerr+"][binIdx[0]:binIdx[-1]+1,:] = temp[:,:,5]

    elif ddict[hname].type()=='Counter':
        f["sumw"][      binIdx[0]:binIdx[-1]+1,:] = temp[:,:,0]
        f["sumw2"][     binIdx[0]:binIdx[-1]+1,:] = temp[:,:,1]
        f["numEntries"][binIdx[0]:binIdx[-1]+1,:] = temp[:,:,2]
    else:
        raise Exception("yikes")


# def namesFromList(L):
    # return [x for x in L if not "/RAW" in x]

# def namesFromListWithRegEx(L, p):
    # return [x for x in L if not p.match(x)]

def mkBinids(hdict):
    binids= []
    for num, hname in enumerate(sorted(list(hdict.keys()))):
        if hname.endswith("]"): continue
        ao = hdict[hname]
        base = ao.path().split("[")[0].replace("/","|")
        if ao.type()=="Scatter1D" or ao.type()=="Scatter2D":
            temp = ["{}#{}".format(base, i) for i in range(len(ao))]
        elif ao.type()=="Counter":
            temp = ["{}#{}".format(base, 0)]
        else:
            suffixes = ["T", "O", "U"]
            if ao.type() == "Counter":
                suffixes.append(0)
            else:
                suffixes.extend([i for i in range(len(ao))])
            temp = ["{}#{}".format(base, s) for s in suffixes]

        binids.extend(temp)
    return binids



D = yoda.readYODA(sys.argv[1])
L = sorted(list(D.keys()))

import re
# p = re.compile("^\/RAW")

names = [x for x in L ]# if not "/RAW" in x]
central     = [x for x in names if not x.endswith("]")]
variations  = [x for x in names if     x.endswith("]")]

# In principle one probably should check that all variations are always the
# same, we assume this is the case vere here
var = []
for c in central:
    var.append([x for x in variations if x.startswith(c+"[")])

## Thats the weight and weight variation order
VVV = ["CentralWeight"]
p=re.compile("\[(.*?)\]")
for x in var[0]:
    try:
        VVV.append(p.findall(x)[0])
    except:
        print(x)

binids = binids = mkBinids(D)

f = h5py.File("Rivet.h5", "w")
createDatasets(f, binids, VVV )

hname=central[0]
for hname in central:
    fillDatasets(f, binids, VVV, D, hname)

f.close()
