import pykep as _pk
import numpy as _np

def add_mit(ax,
    mit,
    mu,
    units=_pk.AU,
    N=60,
    c_segments=["royalblue", "indianred"],
    figsize=(5, 5),
    **kwargs
):
    """
    Plot a Multiple Impulse Trajectory (mit) stored in the mit format:

    mit = [ [[r,v], DV, DT], ... ]


    Args:
        *mit* (:class:`list`): The decision vector in the correct tof encoding.

        *mu* (:class:`float`): The gravitational parameter

        *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`, optional): The 3D axis to plot on. Defaults to None.

        *units* (:class:`float`, optional): The unit scale for the plot. Defaults to pk.AU.

        *N* (:class:`int`, optional): The number of points to use when plotting the trajectory. Defaults to 60.

        *c_segments* (:class:`list`, optional): The colors to alternate the various trajectory segments (inbetween DSMs). Defaults to ["royalblue", "indianred"].

        *figsize* (:class:`tuple`): The figure size (only used if *ax* is None and axis have to be created.), Defaults to (5, 5).

        *\\*\\*kwargs*: Additional keyword arguments to pass to the trajectory plot (common to Lambert arcs and ballistic arcs)

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: The 3D axis where the trajectory was plotted.
    """
    if ax is None:
        ax = _pk.plot.make_3Daxis(figsize=figsize)

    DVs = [_np.linalg.norm(node[1]) for node in mit]
    maxDV = max(DVs)
    DVs = [s / maxDV * 30 for s in DVs]

    # 3 - We loop across grid nodes
    for i, node in enumerate(mit):
        ax.scatter(
            node[0][0][0] / units,
            node[0][0][1] / units,
            node[0][0][2] / units,
            color="k",
            s=DVs[i],
        )

        r_after_dsm = node[0][0]
        v_after_dsm = [a + b for a, b in zip(node[0][1], node[1])]
        _pk.plot.add_ballistic_arc(
            ax,
            [r_after_dsm, v_after_dsm],
            node[2],
            mu,
            N=N,
            units=units,
            c=c_segments[i % len(c_segments)],
            **kwargs
        )
    return ax
