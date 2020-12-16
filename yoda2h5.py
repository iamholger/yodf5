#!/usr/bin/env python3


import yoda, sys
import h5py
import numpy as np

def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    # Fix size, sometimes there is spillover
    # TODO: replace with while if problem persists
    if len(out) > num:
        out[-2].extend(out[-1])
        out = out[0:-1]

    if len(out) != num:
        raise Exception("something went wrong in chunkIt, the target size differs from the actual size")

    return out

def createDatasets(f, binids, variations, depth=1, compression=4):
    """
    Create data sets in the HDF5 file.
    """
    nbins=len(binids)
    nvars=len(variations)

    # The fundamental moments/elements of yoda objecs
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

    # The datasets have 3 axes: binid, weight variation, point in parameter space
    for df in floats: f.create_dataset(df, (nbins,nvars,depth), maxshape=(None,None,None), dtype='f' , chunks=True, compression=compression)

    # Lookups --- helps when reading data and reconstucting YODA objects
    f.create_group("Histo1D")
    f.create_group("Histo2D")
    f.create_group("Profile1D")
    f.create_group("Counter")
    f.create_group("Scatter1D")
    f.create_group("Scatter2D")

    # This is the one that works well with hdf5 when reading std::string in C++
    dt = h5py.special_dtype(vlen=str)
    # We use these simple lists as lookup tables to associate the elements of the datasets ^^^ with
    # the actual YODA Analysis objects
    import numpy as np
    f.create_dataset("binids",     data=np.array(binids,     dtype=dt))
    f.create_dataset("variations", data=np.array(variations, dtype=dt))


def dbn0ToArray(dbn):
    return np.array([dbn.sumW(), dbn.sumW2(), dbn.numEntries()])

def dbn1ToArray(dbn):
    """
    The try except block deals with the underflow things not having xmin, xmax
    """
    try:
        return np.array([dbn.sumW(), dbn.sumW2(), dbn.sumWX(), dbn.sumWX2(), dbn.numEntries(), dbn.xMin(), dbn.xMax()])
    except:
        return np.array([dbn.sumW(), dbn.sumW2(), dbn.sumWX(), dbn.sumWX2(), dbn.numEntries(), 0, 0])

def H2dbn2ToArray(dbn):
    """
    The try except block deals with the underflow things not having xmin, xmax
    """
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

def mkSafeHname(hname):
    return hname.replace("/","|")

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

def mkIndexDict(datadict, allbinids):
    ret = {'Histo1D':{}, 'Histo2D':{}, 'Profile1D':{}, 'Scatter1D':{}, 'Scatter2D':{}, 'Counter':{}}

    for hname, v in datadict.items():
        _hname=mkSafeHname(hname)
        try:
            ret[datadict[hname].type()][_hname] = [num for num, binid in enumerate(allbinids) if binid.startswith("{}#".format(_hname))]
        except Exception as e:
            print("oops: ", e)

    return ret

def createIndexDS(f, d_idx):
    for dtype, objects in d_idx.items():
        for _hname, binIdx in objects.items():
            f.create_dataset("{}/{}".format(dtype, _hname), data=binIdx , chunks=True)

def fillDatasets(f, binIdx, variations, ddict, hname, depth=0):

    if len(binIdx) ==0:
        print("Warning, no matching binid for {} --- is this one of the raw ratios maybe???".format(hname))
        return

    if ddict[hname].type()=='Histo1D':
        nFields=7
        fdbn = dbn1ToArray
    elif ddict[hname].type()=='Histo2D':
        nFields=12
        fdbn = H2dbn2ToArray
    elif ddict[hname].type()=='Profile1D':
        fdbn = dbn2ToArray
        nFields=9
    elif ddict[hname].type()=='Scatter2D':
        fdbn = point2DToArray
        nFields=6
    elif ddict[hname].type()=='Scatter1D':
        fdbn = point1DToArray
        nFields=3
    elif ddict[hname].type()=='Counter':
        nFields=3
    else:
        raise Exception("type {} Not implemented".format(ddict[hname].type()))

    # Empty array to be filled and written to datasets
    temp = np.zeros((len(binIdx), len(variations), nFields))

    hids = [hname]
    for v in variations[1:]:
        hids.append("{}[{}]".format(hname, v))

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
        f["sumw"][      binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,0]
        f["sumw2"][     binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,1]
        f["sumwx"][     binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,2]
        f["sumwx2"][    binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,3]
        f["numEntries"][binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,4]
        f["xmin"][      binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,5]
        f["xmax"][      binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,6]

    # elif ddict[hname].type()=='Histo2D':
        # f["sumw"][      binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,0]
        # f["sumw2"][     binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,1]
        # f["sumwx"][     binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,2]
        # f["sumwx2"][    binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,3]
        # f["sumwy"][     binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,4]
        # f["sumwy2"][    binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,5]
        # f["sumwxy"][    binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,6]
        # f["numEntries"][binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,7]
        # f["xmin"][      binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,8]
        # f["xmax"][      binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,9]
        # f["ymin"][      binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,10]
        # f["ymax"][      binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,11]

    # elif ddict[hname].type()=='Profile1D':
        # f["sumw"][      binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,0]
        # f["sumw2"][     binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,1]
        # f["sumwx"][     binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,2]
        # f["sumwx2"][    binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,3]
        # f["sumwy"][     binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,4]
        # f["sumwy2"][    binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,5]
        # f["numEntries"][binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,6]
        # f["xmin"][      binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,7]
        # f["xmax"][      binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,8]

    # elif ddict[hname].type()=='Scatter1D':
        # f["xval"][ binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,0]
        # f["xerr-"][binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,1]
        # f["xerr+"][binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,2]

    # elif ddict[hname].type()=='Scatter2D':
        # f["xval"][ binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,0]
        # f["xerr-"][binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,1]
        # f["xerr+"][binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,2]
        # f["yval"][ binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,3]
        # f["yerr-"][binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,4]
        # f["yerr+"][binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,5]

    # elif ddict[hname].type()=='Counter':
        # f["sumw"][      binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,0]
        # f["sumw2"][     binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,1]
        # f["numEntries"][binIdx[0]:binIdx[-1]+1,:,depth] = temp[:,:,2]
    # else:
        # raise Exception("yikes")



if __name__=="__main__":
    import sys
    import optparse, os, sys
    op = optparse.OptionParser(usage=__doc__)
    op.add_option("-v", "--debug", dest="DEBUG", action="store_true", default=False, help="Turn on some debug messages")
    op.add_option("-o", dest="OUTPUT", default="analysisobjects.h5", help="Output HDF5 file (default: %default)")
    opts, args = op.parse_args()

    YODAFILES = args

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    binids, VVV, aix, aix_flat, central = None, None, None, None, None

    if rank==0:
        # TODO if len(args)==1 and os.path.isdir(args[0]) --- hierarchical reading with pnames finding etc

        # Let's assume they are all consistent TODO add robustness
        DATA0 = yoda.readYODA(args[0])
        L = sorted(list(DATA0.keys()))

        names       = [x for x in L ]# if not "/RAW" in x]
        central     = [x for x in names if not x.endswith("]")]
        variations  = [x for x in names if     x.endswith("]")]

        # TODO In principle one probably should check that all variations are always the
        # same, we assume this is the case here
        var = []
        for c in central:
            var.append([x for x in variations if x.startswith(c+"[")])

        ## Thats the weight and weight variation order we store the data in
        VVV = ["CentralWeight"]
        import re
        p=re.compile("\[(.*?)\]")
        for x in var[0]:
            try:
                VVV.append(p.findall(x)[0])
            except Exception as e:
                print(x, e)

        binids = mkBinids(DATA0)

        # Hierarchical, i.e. top layer is the AnalysisObject type
        aix = mkIndexDict(DATA0, binids)

        # Object name as keys and lists of indices as values
        aix_flat = {}
        for k, v in aix.items(): aix_flat.update(v)

    binids   = comm.bcast(binids,   root=0)
    VVV      = comm.bcast(VVV,      root=0)
    aix      = comm.bcast(aix,      root=0)
    aix_flat = comm.bcast(aix_flat, root=0)
    central  = comm.bcast(central,  root=0)

    # NOTE dataset operations are collective
    #  This require h5py to use and H5 that is build with MPI
    try:
        f = h5py.File(opts.OUTPUT, "w", driver='mpio', comm=MPI.COMM_WORLD)
    except:
        f = h5py.File(opts.OUTPUT, "w")
    createDatasets(f, binids, VVV, depth=len(YODAFILES))
    createIndexDS(f, aix)

    rankwork = chunkIt([i for i in range(len(YODAFILES))], size) if rank==0 else None
    rankwork = comm.scatter(rankwork, root=0)

    # This part is MPI trivial
    for num, findex in enumerate(rankwork):
        DATA = yoda.readYODA(YODAFILES[findex])
        for hname in central:
            _hname=mkSafeHname(hname)
            fillDatasets(f, aix_flat[_hname], VVV, DATA, hname, depth=findex)
        if rank==0:
            print("[{}] --- {}/{} complete".format(rank, num, len(rankwork)))
        sys.stdout.flush()

    f.close()
