from scipy.interpolate import interp2d, interp1d
from math import sqrt

_vinfs_A5 = [0., 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 5.75, 6]
_decls_A5 = [-90, -40, -30, -29, -28.5, -20, -10, 0, 10, 20, 28.5, 29, 30, 40, 90]
_data_A5 = [
    [1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1],
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0],
    [1160, 1160, 1100, 1010, 930, 830, 740, 630, 590, 550],
    [2335.0, 2335.0	, 2195.0, 2035.0, 1865.0, 1675.0, 1480.0, 1275.0, 1175.0, 1075.0],
    [2335.0, 2335.0	, 2195.0, 2035.0, 1865.0, 1675.0, 1480.0, 1275.0, 1175.0, 1075.0],
    [2335.0, 2335.0	, 2195.0, 2035.0, 1865.0, 1675.0, 1480.0, 1275.0, 1175.0, 1075.0],
    [2335.0, 2335.0	, 2195.0, 2035.0, 1865.0, 1675.0, 1480.0, 1275.0, 1175.0, 1075.0],
    [2335.0, 2335.0	, 2195.0, 2035.0, 1865.0, 1675.0, 1480.0, 1275.0, 1175.0, 1075.0],
    [2335.0, 2335.0	, 2195.0, 2035.0, 1865.0, 1675.0, 1480.0, 1275.0, 1175.0, 1075.0],
    [2335.0, 2335.0	, 2195.0, 2035.0, 1865.0, 1675.0, 1480.0, 1275.0, 1175.0, 1075.0],
    [1160, 1160, 1100, 1010, 930, 830, 740, 630, 590, 550],
    [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0],
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1]
]

_vinfs_A551 = [sqrt(elem) for elem in [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 40, 45, 50, 55, 60]]
_data_A551 = [5995, 5780, 5570, 5360, 5160, 4965, 4775, 4585, 4405, 4230, 4055, 3890, 3730, 3570, 3420, 3270, 3130, 2995, 2860, 2670, 2380, 2120, 1900, 1695]

_vinfs_SF = [0, 1, 2, 3, 4, 5]
_decls_SF = [-90, -65, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 65, 90]
_data_SF = [
    [1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3],
    [100., 100.00000, 100.00000, 100.00000, 100.00000, 100.00000],
    [1830.50000, 1830.50000, 1815.90000, 1737.70000, 1588.00000, 1344.30000],
    [1910.80000, 1910.80000, 1901.90000, 1819.00000, 1636.40000, 1369.30000],
    [2001.80000, 2001.80000, 1995.30000, 1891.30000, 1673.90000, 1391.90000],
    [2108.80000, 2108.80000, 2088.60000, 1947.90000, 1708.00000, 1409.50000],
    [2204.00000, 2204.00000, 2167.30000, 1995.50000, 1734.50000, 1419.60000],
    [2270.80000, 2270.80000, 2205.80000, 2013.60000, 1745.10000, 1435.20000],
    [2204.70000, 2204.70000, 2133.60000, 1965.40000, 1712.80000, 1413.60000],
    [2087.90000, 2087.90000, 2060.60000, 1917.70000, 1681.10000, 1392.50000],
    [1979.17000, 1979.17000, 1975.40000, 1866.50000, 1649.00000, 1371.70000],
    [1886.90000, 1886.90000, 1882.20000, 1801.00000, 1614.60000, 1350.50000],
    [1805.90000, 1805.90000, 1796.00000, 1722.70000, 1571.60000, 1327.60000],
    [100.00000, 100.00000, 100.00000, 100.00000, 100.00000, 100.00000],
    [1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3]
]

# Ariane 5: data provided to ESOC by Arianespace when they were still considering Ariane launch for ExoMars. 
# Negative mass values have been substituted with exp(m/1000) to avoid problems.
# Declination has been extended up to 90 degrees.
_vinfs_Ariane5 = [0.5, 1. , 1.5, 2. , 2.5, 3. , 3.5, 4. , 4.5, 5. , 5.5, 6. ]
_decls_Ariane5 = [-90, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 90]
_data_Ariane5 = [
    [4.97870684e-02, 3.01973834e-02, 1.83156389e-02, 1.11089965e-02, 6.73794700e-03, 2.47875218e-03, 9.11881966e-04, 3.35462628e-04, 1.23409804e-04, 4.53999298e-05, 1.67017008e-05, 6.14421235e-06],
    [4.97870684e-02, 3.01973834e-02, 1.83156389e-02, 1.11089965e-02, 6.73794700e-03, 2.47875218e-03, 9.11881966e-04, 3.35462628e-04, 1.23409804e-04, 4.53999298e-05, 1.67017008e-05, 6.14421235e-06],
    [8.20849986e-02, 4.97870684e-02, 3.01973834e-02, 1.83156389e-02, 1.11089965e-02, 6.73794700e-03, 2.47875218e-03, 9.11881966e-04, 3.35462628e-04, 1.23409804e-04, 4.53999298e-05, 1.67017008e-05],
    [1.35335283e-01, 8.20849986e-02, 4.97870684e-02, 2.55785799e-02, 1.83156389e-02, 1.11089965e-02, 6.73794700e-03, 2.47875218e-03, 9.11881966e-04, 3.35462628e-04, 1.23409804e-04, 4.53999298e-05],
    [2.23130160e-01, 1.35335283e-01, 1.35335283e-01, 8.20849986e-02, 4.97870684e-02, 1.83156389e-02, 1.11089965e-02, 4.08677144e-03, 2.47875218e-03, 9.11881966e-04, 3.35462628e-04, 1.23409804e-04],
    [3.67879441e-01, 2.63685018e-01, 2.23130160e-01, 1.89001562e-01, 1.35335283e-01, 4.97870684e-02, 1.83156389e-02, 6.73794700e-03, 4.08677144e-03, 2.47875218e-03, 9.11881966e-04, 3.35462628e-04],
    [6.06530660e-01, 5.13759511e-01, 3.67879441e-01, 2.63685018e-01, 1.89001562e-01, 1.35335283e-01, 4.97870684e-02, 1.83156389e-02, 6.73794700e-03, 4.08677144e-03, 2.47875218e-03, 9.11881966e-04],
    [5.00600000e+03, 4.66700000e+03, 6.06530660e-01, 3.67879441e-01, 2.63685018e-01, 1.89001562e-01, 1.35335283e-01, 4.97870684e-02, 1.83156389e-02, 1.11089965e-02, 6.73794700e-03, 2.47875218e-03],
    [5.47400000e+03, 5.19500000e+03, 4.80500000e+03, 4.31600000e+03, 7.16770194e-01, 5.13759511e-01, 3.67879441e-01, 1.35335283e-01, 4.97870684e-02, 1.83156389e-02, 1.11089965e-02, 6.73794700e-03],
    [5.83500000e+03, 5.61500000e+03, 5.29100000e+03, 4.87000000e+03, 4.35900000e+03, 3.77400000e+03, 3.13600000e+03, 3.67879441e-01, 1.35335283e-01, 4.97870684e-02, 1.83156389e-02, 1.83156389e-02],
    [6.07800000e+03, 5.91000000e+03, 5.64800000e+03, 5.29500000e+03, 4.85600000e+03, 4.34000000e+03, 3.76300000e+03, 3.14100000e+03, 3.67879441e-01, 2.23130160e-01, 1.35335283e-01, 4.97870684e-02],
    [6.19100000e+03, 6.05900000e+03, 5.84400000e+03, 5.54900000e+03, 5.18000000e+03, 4.74400000e+03, 4.25100000e+03, 3.71400000e+03, 3.15500000e+03, 6.06530660e-01, 3.67879441e-01, 1.35335283e-01],
    [6.07300000e+03, 5.95300000e+03, 5.74900000e+03, 5.47000000e+03, 5.12700000e+03, 4.72900000e+03, 4.28600000e+03, 3.80800000e+03, 3.30400000e+03, 2.78500000e+03, 2.26000000e+03, 3.67879441e-01],
    [3.67879441e-01, 3.32871084e-01, 3.01194212e-01, 2.72531793e-01, 2.46596964e-01, 2.23130160e-01, 2.01896518e-01, 1.82683524e-01, 1.65298888e-01, 1.49568619e-01, 1.35335283e-01, 1.22456428e-01],
    [1.35335283e-01, 1.22456428e-01, 1.10803158e-01, 1.00258844e-01, 9.07179533e-02, 8.20849986e-02, 7.42735782e-02, 6.72055127e-02, 6.08100626e-02, 5.50232201e-02, 4.97870684e-02, 4.50492024e-02],
    [2.23130160e-01, 2.63685018e-01, 3.67879441e-01, 4.72366553e-01, 6.06530660e-01, 7.78800783e-01, 4.08100000e+03, 3.50900000e+03, 2.89100000e+03, 2.24400000e+03, 3.67879441e-01, 1.35335283e-01],
    [3.67879441e-01, 4.72366553e-01, 6.06530660e-01, 7.78800783e-01, 4.87400000e+03, 4.39100000e+03, 3.83600000e+03, 3.22000000e+03, 2.55900000e+03, 3.67879441e-01, 1.35335283e-01, 4.97870684e-02],
    [6.06530660e-01, 7.78800783e-01, 5.48400000e+03, 5.13400000e+03, 4.69300000e+03, 4.16700000e+03, 3.56300000e+03, 2.89700000e+03, 2.19000000e+03, 2.23130160e-01, 4.97870684e-02, 1.83156389e-02],
    [5.77300000e+03, 5.58900000e+03, 5.30600000e+03, 4.92400000e+03, 4.44300000e+03, 3.86800000e+03, 3.21000000e+03, 3.67879441e-01, 1.35335283e-01, 4.97870684e-02, 1.83156389e-02, 6.73794700e-03],
    [5.65000000e+03, 5.44100000e+03, 5.12400000e+03, 4.69700000e+03, 4.16100000e+03, 3.52200000e+03, 3.67879441e-01, 1.35335283e-01, 4.97870684e-02, 1.83156389e-02, 6.73794700e-03, 2.47875218e-03],
    [5.47700000e+03, 5.23900000e+03, 4.88200000e+03, 4.40100000e+03, 3.79500000e+03, 3.67879441e-01, 1.35335283e-01, 4.97870684e-02, 1.83156389e-02, 6.73794700e-03, 2.47875218e-03, 9.11881966e-04],
    [5.30200000e+03, 5.02100000e+03, 4.60400000e+03, 4.04400000e+03, 3.67879441e-01, 1.35335283e-01, 4.97870684e-02, 1.83156389e-02, 6.73794700e-03, 2.47875218e-03, 9.11881966e-04, 3.35462628e-04],
    [5.30200000e+03, 5.02100000e+03, 4.60400000e+03, 4.04400000e+03, 3.67879441e-01, 1.35335283e-01, 4.97870684e-02, 1.83156389e-02, 6.73794700e-03, 2.47875218e-03, 9.11881966e-04, 3.35462628e-04]
]

class _launchers:
    """
    This class contains a few functions that return the mass launchers can deliver
    to a certain declination / vinf.

    .. note::

       In pykep the object pykep.trajopt.launchers is already an instance of this class and is to be used
       as it has all the data preallocated upon import.

    Examples:
        >>> import pykep as pk
        >>> mass = pk.trajopt.launchers.soyuzf(4.5, 33.21)
        >>> mass2 = pk.trajopt.launchers.atlas501(4.5, 33.21)

    """
    def __init__(self):
        self._atlas501 = interp2d(_vinfs_A5, _decls_A5, _data_A5, kind='linear', fill_value=0.1, copy=False)
        self._atlas551 = interp1d(_vinfs_A551, _data_A551, kind='linear', fill_value=0.1, copy=False, bounds_error=False)
        self._soyuzf = interp2d(_vinfs_SF, _decls_SF, _data_SF, kind='linear', fill_value=1e-3, copy=False)
        self._ariane5 = interp2d(_vinfs_Ariane5, _decls_Ariane5, _data_Ariane5, kind='linear', fill_value=1e-3, copy=False)
    def atlas501(self, vinfs, decls):
        """atlas501(vinfs, decls)

        Computes the mass that the Atlas 501 launcher can deliver to a certain vinf and declination.
        If the inputs are arrays, then a mesh is considered and the mass is returned on points of the mesh

        Args:
            - vinfs (``float`` or array-like): the hyperbolic escape velocity in km/s
            - decls (``float`` or array-like): the declination in degrees

        Returns:
            Numpy array containg the mass delivered to escape with said declinations and magnitudes.
        
        """
        return self._atlas501(vinfs, decls)

    def atlas551(self, vinfs):
        """atlas551(vinfs)

        Computes the mass that the Atlas 551 launcher can deliver to a certain vinf
        If the inputs are arrays, then a mesh is considered and the mass is returned on points of the mesh

        Args:
            - vinfs (``float`` or array-like): the hyperbolic escape velocity in km/s

        Returns:
            Numpy array containg the mass delivered to escape with said magnitudes.

        """
        return self._atlas551(vinfs)

    def soyuzf(self, vinfs, decls):
        """soyuzf(vinfs, decls)

        Computes the mass that the Soyutz-Fregat launcher can deliver to a certain vinf and declination.
        If the inputs are arrays, then a mesh is considered and the mass is returned on points of the mesh

        Args:
            - vinfs (``float`` or array-like): the hyperbolic escape velocity in km/s
            - decls (``float`` or array-like): the declination in degrees

        Returns:
            Numpy array containg the mass delivered to escape with said declinations and magnitudes.
        """
        return self._soyuzf(vinfs, decls)
    def ariane5(self, vinfs, decls):
        """ariane5(vinfs, decls)

        Computes the mass that the Ariane5 launcher can deliver to a certain vinf and declination, assuming
        a launch from Kourou. Data provided to ESOC by Arianespace when Ariane launch for ExoMars was an option.
        If the inputs are arrays, then a mesh is considered and the mass is returned on points of the mesh.

        Args:
            - vinfs (``float`` or array-like): the hyperbolic escape velocity in km/s
            - decls (``float`` or array-like): the declination in degrees

        Returns:
            Numpy array containg the mass delivered to escape with said declinations and magnitudes.
        """
        return self._ariane5(vinfs, decls)
