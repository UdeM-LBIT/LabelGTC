#ifndef GENESUBTREECORRECTOR_H
#define GENESUBTREECORRECTOR_H

#include "trees/node.h"
#include "trees/genespeciestreeutil.h"
#include "supergenetreemaker.h"


#include <unordered_map>
#include <vector>
#include <algorithm>
#include <set>

using namespace std;

class GeneSubtreeCorrector
{
public:
    GeneSubtreeCorrector();

    Node* GetSubtreeTripletRespectingHistory(Node* geneTree, Node* speciesTree,
                                             unordered_map<Node *, Node *> &lca_mapping,
                                             vector<Node *> &geneSubtreeRoots, bool mustRetainDupSpec);


private:

    Node* minTRS_Rec(Node* geneTree, Node* speciesTree, unordered_map<Node *, Node *> &lca_mapping,
                     set<Node*> &markedGeneNodes, bool mustRetainDupSpec);

};

#endif // GENESUBTREECORRECTOR_H
