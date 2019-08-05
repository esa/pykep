class dbscan():
    """
    This class can be used to locate areas of the interplanetsry space that are 'dense' at one epoch.
    Essentially, it locates planet clusters
    """
    from pykep.core import AU, EARTH_VELOCITY

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

        - planet_list = list of pykep planets (typically thousands)
        """
        self._asteroids = planet_list
        self.labels = None
        self.n_clusters = None
        self.members = None
        self.core_members = None

    def _orbital_metric(self, r, v):
        from pykep.core import DAY2SEC
        DV2 = [a / (self._T * DAY2SEC) for a in r]
        DV1 = [a + b for a, b in zip(DV2, v)]
        return DV1 + DV2

    def cluster(self, t, eps=0.125, min_samples=10, metric='orbital', T=180, ref_r=AU, ref_v=EARTH_VELOCITY):
        """
        USAGE: cl.cluster(t, eps=0.125, min_samples=10, metric='orbital', T=180, ref_r=AU, ref_v=EARTH_VELOCITY):

        - t: epoch (in MJD2000)
        - eps: max distance between points in a cluster
        - min_samples: minimum number of samples per cluster
        - metric: one of 'euclidean', 'euclidean_r', orbital'
        - T: average transfer time (used in the definition of the 'orbital' metric)
        - ref_r         reference radius   (used as a scaling factor for r if the metric is 'euclidean' or 'euclidean_r')
        - ref_v         reference velocity (used as a scaling factor for v if the metric is 'euclidean')
        """
        import pykep
        import numpy
        from sklearn.cluster import DBSCAN

        self._epoch = pykep.epoch(t)

        if metric == 'euclidean':
            self._X = [
                [elem for tupl in p.eph(self._epoch) for elem in tupl] for p in self._asteroids]
            scaling_vector = [ref_r] * 3
            scaling_vector += [ref_v] * 3
        elif metric == 'euclidean_r':
            self._X = [list(p.eph(self._epoch)[0]) for p in self._asteroids]
            scaling_vector = [ref_r] * 3
        elif metric == 'orbital':
            self._T = T
            self._X = [self._orbital_metric(
                *p.eph(self._epoch)) for p in self._asteroids]
            scaling_vector = [1.] * 6  # no scaling
        self._X = numpy.array(self._X)

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
            from pykep.orbit_plots import plot_planet
            members = self.core_members if only_core else self.members
            for label in members if clusters is None else clusters:
                for planet in members[label]:
                    plot_planet(
                        self._asteroids[planet], t0=self._epoch, s=0, axes=axis)

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
            print(
                "cluster_id should be larger then 0 and smaller than the number of clusters (-1)")
            return
        if len(epochs) != 9:
            print("The epochs requested must be exactly 9 as to assemble 3x3 subplots")
            return

        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        from pykep.orbit_plots import plot_planet
        from pykep import epoch

        if only_core:
            ids = self.core_members[cluster_id]
        else:
            ids = self.members[cluster_id]

        fig = plt.figure()
        for i, ep in enumerate(epochs):
            axis = fig.add_subplot(3, 3, i + 1, projection='3d')

            plt.axis('off')
            plt.title(epoch(ep).__repr__()[:11])
            for pl in self._asteroids[::skip]:
                axis = plot_planet(pl, axes=axis, alpha=0.05, s=0)
            for cluster_member in ids:
                r, _ = self._asteroids[cluster_member].eph(epoch(ep))
                axis.scatter([r[0]], [r[1]], [r[2]], marker='o', alpha=alpha)

        plt.draw()
        plt.show()
        return fig
