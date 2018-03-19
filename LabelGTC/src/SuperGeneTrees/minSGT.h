
#include "trees/node.h"
#include "trees/newicklex.h"
#include "trees/genespeciestreeutil.h"
#include "supergenetreemaker.h"
#include "genesubtreecorrector.h"
#include "trees/polysolver.h"

#include <iostream>
#include <string>
#include <map>
#include <time.h>
#include "div/tinydir.h"

using namespace std;

string DoSuperGeneTree(string gcontent, string scontent, bool preserveDupSpec, string clades_to_preserve, string treated_trees, string outputmode="", int limit=1);
string DoSubtreeCorrection(string gcontent, string scontent, bool preserveDupSpec, string markedNodesMode = "", string outputmode = "tree");

string DoSGTOnHighSpecs(string gcontent, string scontent, string clades_to_preserve, string treated_trees, bool preserveDupSpec, int maxNBLeaves = 999999, int maxNBSubtrees = 999999);
void DoPolytomyCorrection(string gcontent, string scontent);

/*int minSGT(int argc, char *argv[]);
{

    here I store my previous test command line args
    -p -gn "((A1__A, C1__C),(B1__B, (C2__C,D2__D)));((B2__B,D3__D),(B1__B,C3__C));" -sn "((A,B),(C,D))"
    -p -m sub -gn "(((A1__A, B1__B)m, (A2__A, B2__B)m), ((C3__C, D3__D), ((A4__A,B4__B),(C4__C,D4__D)))m);" -sn "((A,B),(C,D))"

    -m sub -sn "(((A,B),C), (D,E));" -gn "(((A1__A, B1__B)m,C1__Cm), ((D1__D,E1__E),(((A2__A,B2__B),C2__C), (D2__D,E2__E)))m);"

    -m sub -sn "(((A,B),C), (D,E));" -gn "(((((A1__A, D1__D),B1__B)m,((B2__B, C2__C),E2__E)m),(((D3__D, E3__E),(D4__D, E4__E))m,((A5__A, B5__B), (C5__C, C6__C))m)),((A6__A, D6__D),(B7__B, (B8__B, C8__C)))m);"

    BATCH
    -m sub -s "C:/Users/Manuel/Desktop/projects/code/data/species_tree.ensembl.topology.nw" -b "C:/Users/Manuel/Desktop/projects/code/data/gene_trees_may2016/"

    ENSGT00460000041710.txt
    -m sub -p -s "C:/Users/Manuel/Desktop/projects/code/data/species_tree.ensembl.topology.nw" -g "C:/Users/Manuel/Desktop/projects/code/data/gene_trees_may2016/ENSGT00460000041710.txt"


    FUCK THIS BELOW
    -m sub -s "C:/Users/Manuel/Desktop/projects/code/data/species_tree.ensembl.topology.nw" -b "C:/Users/Manuel/Desktop/projects/code/data/test_gene_trees/"

    -m sub -s "C:/Users/Manuel/Desktop/projects/code/data/species_tree.ensembl.topology.nw" -b "C:/Users/Manuel/Desktop/projects/code/data/test_gene_trees/"

    //gene tree label format is GENENAME__SPECIESNAME

    //below = supertree mode

}*/
