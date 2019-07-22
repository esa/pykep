import numpy as np


def load_gravity_model(fp):
    """Load spherical harmonics model.

    Args:
        - fp (``str``): path to model file.

    Returns:
        - planet_radius (``float``): Radius of central body (km)
        - gravity_param (``float``): Gravitational parameter of central body (km3/s2)
        - c_coefficients (``array-like``): Array of C_(n, m) coefficients.
        - s_coefficients (``array-like``): Array of S_(n, m) coefficients.
        - degree (``int``): Maximum degree of model.
        - order (``int``): Maximum order of model.

    Example::

        r, mu, c, s, n, m = pykep.util.load_gravity_model('gravity_models/Earth/egm96.txt')

    File format:
        The input file must follow a specific format. The first row contains planetary radius in km and
        gravitational parameter in km3/s2. The second row until the end contains normalised C and S 
        coefficients for every degree and order. The file is delimited using ';'.
    
    File::

          1| radius;mu;degree;order
          2| 0;0;C(0,0);S(0,0)
          3| 1;0;C(1,0);S(1,0)
          4| ...
    """
    # todo: check file

    planet_radius, gravity_param, degree, order, array = _read_model_file(fp)

    c_coefficients, s_coefficients = _convert_data(array, degree, order)

    return planet_radius, gravity_param, c_coefficients, s_coefficients, degree, order


def _read_model_file(fp):
    """Read model data from file"""
    with open(fp) as f:
        header = np.genfromtxt(f, delimiter=";", max_rows=1)

        if len(header) != 4:
            raise IndexError("file header needs at least 4 values: radius, mu, degree, order. " +
                             f"Header contains ({len(header)}) values.")

        radius = float(header[0])
        gravity_const = float(header[1])
        degree = int(header[2])
        order = int(header[3])

        array = np.genfromtxt(f, delimiter=";", skip_header=1)

        return radius, gravity_const, degree, order, array


def _convert_data(array, max_degree, max_order):
    """Convert array from file to two arrays containing only C or S coefficients."""
    c_array = np.zeros((max_degree + 1, max_order + 1))
    s_array = np.zeros((max_degree + 1, max_order + 1))

    for row in array:
        degree = int(row[0])
        order = int(row[1])

        c = row[2]
        s = row[3]

        c_array[degree, order] = c
        s_array[degree, order] = s

    return c_array, s_array


if __name__ == "__main__":
    load_coefficients("C:/Users/Bert van den Abbeele/Documents/Python/HarmonicsPropagator/tests/jgl150q1.txt")
