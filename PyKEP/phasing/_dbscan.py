import PyKEP
import numpy


class dbscan():

    """
    This class can be used to locate areas of the interplanetsry space that are 'dense' at one epoch.
    Essentially, it locates planet clusters
    """

    def _axis_equal_3d(self, ax):
        """Rescales 3D axis limits using equal scale."""
        import numpy
        extents = numpy.array(
            [getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
        sz = extents[:, 1] - extents[:, 0]
        centers = numpy.mean(extents, axis=1)
        maxsize = max(abs(sz))
        r = maxsize / 2
        for ctr, dim in zip(centers, 'xyz'):
            getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

    def __init__(self, planet_list):
        """
        USAGE: cl = dbscan(planet_list):

        - planet_list = list of PyKEP planets (typically thousands)
        """
        self._asteroids = planet_list
        self.labels = None
        self.n_clusters = None
        self.members = None
        self.core_members = None

    def cluster(self, t, with_velocity=True, scaling='astro', eps=0.125, min_samples=10):
        """
        USAGE: cl.cluster(self, t, with_velocity=True, scaling='astro', eps=0.125, min_samples=10)

        - t: epoch (in MJD2000)
        - with_velocity: when True clusters by position and velocity, otherwise only position is used
        - scaling: one of
          - None, or
          - 'standard' (removing mean and scale to standard variance), or
          - 'astro' (scaling by PyKEP.AU and PyKEP.EARTH_VELOCITY)
        - eps: max distance between points in a cluster
        - min_samples: minimum number of samples per cluster
        """
        import PyKEP
        import numpy
        from sklearn.preprocessing import StandardScaler
        from sklearn.cluster import DBSCAN

        self._scaling = scaling
        self._epoch = PyKEP.epoch(t)

        if with_velocity:
            self._X = [
                [elem for tupl in p.eph(self._epoch) for elem in tupl] for p in self._asteroids]
        else:
            self._X = [list(p.eph(self._epoch)[0]) for p in self._asteroids]
        self._X = numpy.array(self._X)

        self._scaler = None
        if self._scaling == 'standard':
            self._scaler = StandardScaler().fit(self._X)
            self._X = self._scaler.transform(self._X)
        elif self._scaling == 'astro':
            scaling_vector = [PyKEP.AU] * 3
            if with_velocity:
                scaling_vector += [PyKEP.EARTH_VELOCITY] * 3
            scaling_vector = numpy.array(scaling_vector)
            self._X = self._X / scaling_vector[None, :]

        self._db = DBSCAN(eps=eps, min_samples=min_samples).fit(self._X)
        self._core_samples = self._db.core_sample_indices_

        self.labels = self._db.labels_
        self.n_clusters = len(
            set(self.labels)) - (1 if -1 in self.labels else 0)

        self.members = {}
        self.core_members = {}
        for label in set(self.labels):
            if int(label) == -1:
                continue
            self.members[int(label)] = [index[0]
                                        for index in numpy.argwhere(self.labels == label)]
            self.core_members[int(label)] = [
                index for index in self._core_samples if self.labels[index] == label]

        if self._scaling == 'standard':
            self._X = self._scaler.inverse_transform(self._X)
        elif self._scaling == 'astro':
            self._X = self._X * scaling_vector[None, :]

    def pretty(self):
        """Prints the cluster lists."""
        if self.labels is None:
            return
        print("Number of clusters: %d" % self.n_clusters)
        print("Size of dataset: %s" % str(self._X.shape))
        print("Scaling: %s" % str(self._scaling))
        print("Epoch: %s" % str(self._epoch))
        for label in list(self.members.keys()):
            print("cluster %d (%d - %d): %s" % (label, len(self.members[label]),
                                                len(self.core_members[label]), str(self.members[label])))

    def plot(self, ax=None, clusters=None, orbits=False, only_core=False):
        """Plots the clusters."""
        if self.n_clusters < 1:
            return

        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D

        if ax is None:
            fig = plt.figure()
            axis = fig.add_subplot(111, projection='3d')
        else:
            axis = ax

        axis.view_init(elev=30.0, azim=135.0)
        axis.set_aspect('equal')

        if orbits:
            from PyKEP.orbit_plots import plot_planet
            members = self.core_members if only_core else self.members
            for label in members if clusters is None else clusters:
                for planet in members[label]:
                    plot_planet(
                        self._asteroids[planet], t0=self._epoch, s=0, ax=axis)

        X, labels = list(zip(*[(x, label) for (x, label) in zip(self._X, self.labels)
                               if label > -.5 and (clusters is None or label in clusters)]))
        data = [[x[0], x[1], x[2]] for x in X]
        axis.scatter(*list(zip(*data)), c=labels, alpha=0.5)

        self._axis_equal_3d(axis)

        if ax is None:
            plt.show()
        return axis

    def plot_cluster_evolution(self, cluster_id=None, only_core=False, epochs=range(7500, 8400, 100), skip=100, alpha=0.3):
        """
        Plots a cluster evolution at 9 prefixed epochs.


        """
        if self.n_clusters < 1:
            print("No clusters have been found yet")
            return
        if cluster_id >= self.n_clusters or cluster_id < 0:
            print("cluster_id should be larger then 0 and smaller than the number of clusters (-1)")
            return
        if len(epochs) != 9:
            print("The epochs requested must be exactly 9 as to assemble 3x3 subplots")
            return

        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        from PyKEP.orbit_plots import plot_planet
        from PyKEP import epoch

        if only_core:
            ids = self.core_members[cluster_id]
        else:
            ids = self.members[cluster_id]

        fig = plt.figure()
        for i, ep in enumerate(epochs):
            print i
            axis = fig.add_subplot(3, 3, i + 1, projection='3d')

            plt.axis('off')
            plt.title(epoch(ep).__repr__()[:11])
            for pl in self._asteroids[::skip]:
                axis = plot_planet(pl, ax=axis, alpha=0.05, s=0)
            for cluster_member in ids:
                r, v = self._asteroids[cluster_member].eph(epoch(ep))
                axis.scatter([r[0]], [r[1]], [r[2]], marker='o', alpha=alpha)


        plt.draw()
        plt.show()
        return fig

del PyKEP, numpy
