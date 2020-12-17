import h5py

def read(fname):
    with h5py.File(fname, "r") as f:
        data = f["Histo1D"][:]
    return data

if __name__ == "__main__":
    import sys

    #TODO scaling and stuff

    fout = sys.argv[1]

    tomerge = sys.argv[2:]

    assert(fout not in tomerge)
    import os
    assert(not os.path.exists(fout))

    base = read(tomerge[0])

    for fn in tomerge[1:]:
        base += read(fn)

    with h5py.File(fout, "w") as f:
        f.create_dataset("Histo1D", data=base, compression="gzip", compression_opts=1)

