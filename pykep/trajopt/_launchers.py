from scipy.interpolate import interp2d

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
        self._soyuzf = interp2d(_vinfs_SF, _decls_SF, _data_SF, kind='linear', fill_value=1e-3, copy=False)
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
