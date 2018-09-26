import numpy as np

header_map = {}
header_map["ARMA_MAT_BIN_IS004"] = "<i4"
header_map["ARMA_MAT_BIN_FN008"] = "<f8"
header_map["ARMA_MAT_BIN_FN004"] = "<f4"

header_map_copy = header_map.copy()

for key, value in header_map_copy.iteritems():
    header_map[value] = key

def load_mat(filename, verbose=False):
    f = open(filename, "rb")
    header = f.readline().replace('\r', '').replace('\n', '')
    if verbose:
        print "Reading file: %s"%filename
    try:
        dtype = header_map[header]
        if verbose:
            print "Header: %s"%header
    except:
        raise ValueError("Datatype not understood!")
    else:
        sizes = f.readline()
        n_rows = int(sizes.split()[0])
        n_cols = int(sizes.split()[1])
        if verbose:
            print "Matrix size: %s x %s"%(n_rows, n_cols)
        x = np.fromfile(f, dtype=dtype)
        assert x.size == n_rows*n_cols
        x = np.reshape(x, [n_rows, n_cols], order="F")
    f.close()
    return x

def save_mat(x, filename, verbose=False):
    n, m = x.shape
    if verbose:
        print "Writing array of size: %i x %i"%(n, m)
    try:
        header = header_map[x.dtype.str]
        if verbose:
            print "Header: %s"%header
    except:
        raise ValueError("Datatype not understood!")

    f = open(filename, "wb")
    f.write("%s\n"%header)
    f.write("%s %s\n"%(n, m))
    xx = np.reshape(x, [n*m, 1], order="F")
    xx.tofile(f)
    f.close()

if __name__ == "__main__":
    x = load_mat("test.bin", verbose=False)
    save_mat(x, "test_w.bin", verbose=False)
