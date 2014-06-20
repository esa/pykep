import numpy as np
import PyKEP as pk



def _eph_normalize( eph, ref_r=pk.AU, ref_v=pk.EARTH_VELOCITY ):
	"""
	Normalize a body's ephemerides:
	    ( r / refn_r, v / refn_v )
	where `refn_r` and `refn_v` are reference values you can set in the
	arguments (defaults: refn_r == pk.AU, refn_v == pk.EARTH_VELOCITY).
	
	Accepts either a single ephemeride, or a matrix with one per row.
	
	Example:
	>>> eph_normalize( asteroid_list[0].eph(TRAJ_START_MIN) )
	array([ -1.80619704e-01,   9.66554711e-01,  -1.46653866e-05,
	        -9.99416814e-01,  -1.87484002e-01,   4.12493855e-06])
	"""
	if type(eph) == tuple:
		eph = np.hstack( eph )
	e = eph.reshape(-1,6)
	
	# normalize r
	e[:, :3] /= ref_r
	# normalize v
	e[:,-3:] /= ref_v
	
	return e.reshape( eph.shape )
	
#--------------------------------------# Spatial indexing

def make_kdtree( planets_list, t, ref_r=pk.AU, ref_v=pk.EARTH_VELOCITY ):
	"""
	Returns a kd-tree data structure indexing the normalized ephemerides
	of planets_list, at epoch t. 

	make_kdtree( planets_list, t, ref_r=AU, ref_v=EARTH_VELOCITY )

	- planets_list: list of PyKEP.planet objects
	- t: epoch
	- ref_r: non dimensional units for distance 
	- ref_v: non dimensional units for velocities

	The returned kd-tree can then be used for efficient nearest-neighbor queries.
	
	See also:
		http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html
		https://en.wikipedia.org/wiki/K-d_tree

	Examples::

	  ast_list = [PyKEP.planet_gtoc7(idx) for idx in range(1,16256)]
	  kdt = PyKEP.phasing.make_kdtree(ast_list,ePyKEP.epoch(6754.3))
	"""

	from scipy.spatial import cKDTree
	

	# asteroids list converted to array of object, 
	# for efficient selection of multiple indices (as well as faster indexing)
	planets_list = np.array( planets_list, dtype=np.object )
	
	# get asteroids' ephemerides into a matrix of shape (len(planets_list), 2, 3) 
	e = np.array( [ a.eph( t ) for a in planets_list ] )

	# reshape memory area, so each asteroid's ephemeride gets represented by
	# a single 6 dimensional vector
	e = e.reshape( (e.shape[0], -1) )

	# normalize the full matrix
	# (if the `normalize_f` is set to None, no normalization takes place)
	_eph_normalize( e, ref_r, ref_v )
	
	return cKDTree( e ), planets_list, t
	
	
#--------------------------------------# Spatial querying

def find_neighbours( (kdtree, planets_list, t), query_planet, query_type='knn', ref_r=pk.AU, ref_v=pk.EARTH_VELOCITY, *args, **kwargs ):
	"""
	Finds the neighbours of a given planet at a given epoch. The user may query for the 
	k-nearest neighbours or for all neighbours within a given "distance"

	find_neighbours( kdtree, query_planet, query_type='knn', ref_r=pk.AU, ref_v=pk.EARTH_VELOCITY, *args, **kwargs ):

	- kdtree: a tuple (kdtree, planets_list, t) returned by PyKEP.phasing.make_kdtree
	- query_planet: the planet we want to find neighbours of
	- query_type: one of 'knn' or 'ball'.
	- ref_r: non dimensional units for distance 
	- ref_v: non dimensional units for velocities
	- *args, **args: according to the query type

	Returns (neighb, neighb_ids, dists), where dist is only computed if 'knn' type query is made

	The following kinds of spatial queries are currently implemented:
	
	query_type = 'knn':
		Obtain the `k` nearest asteroids
		For arguments, see:
		http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query.html
	
	query_type = 'ball':
		Obtain all asteroids within a radius `r` of the given query asteroid
		For arguments, see:
		http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query_ball_point.html
	"""
	if type(query_planet) == int:
		query_planet = planets_list[query_planet]
	
	# generate the query vector
	x = query_planet.eph( t )
	x = _eph_normalize( x, ref_r, ref_v )
	
	
	if query_type == 'knn':
		# Query for the k nearest neighbors
		# http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query.html
		dists, idxs = kdtree.query( x, *args, **kwargs )
	elif query_type == 'ball':
		# Query for all neighbors within a sphere of given radius
		# http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query_ball_point.html
		idxs = kdtree.query_ball_point( x, *args, **kwargs )
		dists = [None] * len(idxs)
	else:
		raise Exception( 'Unrecognized query type: %s' % str(query_type) )
	
	
	neighb = [
		( planets_list[i], i, d )		# (ast. object, ast. ID, distance)
		for i,d in zip(idxs, dists)
		]
		
	# split into three lists, one of objects, one of IDs, and one for distances
	neighb, neighb_ids, dists = zip( *neighb ) if neighb != [] else ([], [], [])
	
	return neighb, neighb_ids, dists

del pk, np
	