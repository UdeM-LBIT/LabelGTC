#ifndef SUPERGENETREEMAKER_H
#define SUPERGENETREEMAKER_H

#include <string>
#include <vector>
#include <unordered_map>

#include "trees/node.h"
#include "div/util.h"
#include "trees/genespeciestreeutil.h"


using namespace std;

//NOTE: trees must have internal nodes labeled uniquely to use this class
class TreeLabelIntersectionInfo
{
public:
    TreeLabelIntersectionInfo();

    void ComputeAllIntersections(vector<Node*> trees);

    void ComputeIntersections(Node *tree1, Node *tree2);

    void AddIntersection(string lbl1, string lbl2);

    bool Intersect(string lbl1, string lbl2);

    bool Intersect(vector<Node*> trees);

    bool IsPartitionIntersecting(vector<Node*> treesLeft, vector<Node*> treesRight);

private:
    unordered_set<string> intersecting_nodes;   //key is inter1-inter2 (separated by this)
    string key_separator;



};

class SuperGeneTreeMaker
{
public:
    SuperGeneTreeMaker();

    //returns a supertree + DL cost
    // En changed this to a vector
    vector<pair<Node*, int>> GetSuperGeneTreeMinDL(vector<Node *> &trees, vector<Node *> &clades_to_preserve, vector<Node *> &treated_trees, vector<unordered_map<Node *, Node *> > &lca_mappings,
                                       Node* speciesTree, bool mustPreserveDupSpec, bool isFirstCall = true, int limit = 1);

private:


    void ApplyNextConfig(vector<int> &counters);

    TreeLabelIntersectionInfo intersectionInfo;
    // En changed this to a map with a list instead of a pair
    unordered_map<string, vector<pair<Node*, int>> > recursionCache;


};

#endif // SUPERGENETREEMAKER_H
