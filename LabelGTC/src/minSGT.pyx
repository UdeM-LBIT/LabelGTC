from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "SuperGeneTrees/minSGT.h":
	string DoSuperGeneTree(string gcontent, string scontent, bool preserveDupSpec, string outputmode, string clades_to_preserve, string treated_trees)


cpdef getMinSGT(string gcontent, string scontent, bool preserveDupSpec, string clades, string trees, string outmode=""):


	res = DoSuperGeneTree(gcontent, scontent, preserveDupSpec, clades, trees, outmode)
	return res
