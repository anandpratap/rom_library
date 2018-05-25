cimport rom

cdef class gemsrom1:

    cdef rom.gemsrom *_c_gemsrom
    def __cinit__(self):
        self._c_gemsrom = rom.create_gemsrom()

    def initialize(self):
        rom.gemsrom_initialize(self._c_gemsrom, 0)

cdef class mainrom1:
    cdef rom.mainrom *_c_mainrom
    def __cinit__(self):
        self._c_mainrom = rom.create_mainrom()

    def calc_svd(self):
        rom.mainrom_calc_svd(self._c_mainrom)
