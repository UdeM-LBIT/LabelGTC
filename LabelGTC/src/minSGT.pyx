from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "SuperGeneTrees/minSGT.h":
	string DoSuperGeneTree(string gcontent, string scontent, bool preserveDupSpec, string clades_to_preserve, string treated_trees, string outputmode, int limit)


cpdef getMinSGT(string gcontent, string scontent, bool preserveDupSpec, string clades, string trees, string outmode="", int limit=1):

	res = DoSuperGeneTree(gcontent, scontent, preserveDupSpec, clades, trees, outmode, limit)
	return res