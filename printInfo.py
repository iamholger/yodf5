import h5py

def readAttribute(fname, dsname, key):
    with h5py.File(sys.argv[1], "r") as f:
        ddd = dict(f[dsname].attrs)
    return ddd[key]


if __name__ == "__main__":
    import sys

    with h5py.File(sys.argv[1], "r") as f:
        for aot, ds in f.items():
            print("Got AnalysisObject type {}".format(aot))
            print("Dataset is of shape {}".format(ds.shape))
            print("Attributes:")
            for k, v in ds.attrs.items():
                print("Atribute {} with len {} ".format(k, len(v)))

    print(readAttribute(sys.argv[1], "Histo1D", "binids"))
