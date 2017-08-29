#include "genesubtreecorrector.h"

GeneSubtreeCorrector::GeneSubtreeCorrector()
{

}


Node* GeneSubtreeCorrector::minTRS_Rec(Node* geneTree, Node* speciesTree, unordered_map<Node *, Node *> &lca_mapping,
                 set<Node*> &markedGeneNodes, bool mustRetainDupSpec)
{

    if (geneTree->IsLeaf() || markedGeneNodes.find(geneTree) != markedGeneNodes.end())
    {
        Node* ret = new Node(false);
        ret->CopyFrom(geneTree);
        lca_mapping[ret] = lca_mapping[geneTree];
        return geneTree;
    }

    Node* ch1 = geneTree->GetChild(0);
    Node* ch2 = geneTree->GetChild(1);

    //both children are marked === here is a cherry
    //so we call minSGT on that
    if (markedGeneNodes.find(ch1) != markedGeneNodes.end()
            && markedGeneNodes.find(ch2) != markedGeneNodes.end() )
    {
        vector<Node*> subtrees;
        subtrees.push_back(ch1);
        subtrees.push_back(ch2);
        vector< unordered_map<Node*, Node*> > subtrees_lca_mapping;

        //why twice the same line below?  because supergenetreemaker needs one lca mapping
        //per gene tree.  But in our case the same lca mapping can be used for both
        subtrees_lca_mapping.push_back(lca_mapping);
        subtrees_lca_mapping.push_back(lca_mapping);


        SuperGeneTreeMaker supermaker;

        vector<Node *> clades_to_preserve = {};
        vector<Node *> treated_trees = {};

        pair<Node*, int> res = supermaker.GetSuperGeneTreeMinDL(subtrees, clades_to_preserve, treated_trees, subtrees_lca_mapping, speciesTree, mustRetainDupSpec);


        Node* newsubtree = res.first;

        lca_mapping[newsubtree] = lca_mapping[geneTree];


        return newsubtree;

    }
    //ch1 is marked, not ch2 (or converse) ==> try grafting minRTS(ch2) on every branch of ch1
    else if (markedGeneNodes.find(ch1) != markedGeneNodes.end() ||
             markedGeneNodes.find(ch2) != markedGeneNodes.end())
    {

        Node* markedTree = ch1;
        Node* otherTree = ch2;
        if (markedGeneNodes.find(ch2) != markedGeneNodes.end())
        {
            markedTree = ch2;
            otherTree = ch1;
        }

        //we will work with copies
        //here, because we're in a hurry to submit, we just try every possibility,
        //compute the reconciliation cost on each possibility
        Node* markedTree2 = new Node(false);
        markedTree2->CopyFrom(markedTree);

        Node* otherTree2 = minTRS_Rec(otherTree, speciesTree, lca_mapping, markedGeneNodes, mustRetainDupSpec);


        int minDL = 99999;
        Node* minNode = NULL;

        //attempt 1: original topology
        Node* tmp = new Node(false);
        tmp->AddSubTree(markedTree2);
        tmp->AddSubTree(otherTree2);

        unordered_map<Node*, Node*> lca_mapping_tmp =
                GeneSpeciesTreeUtil::Instance()->GetLCAMapping(tmp, speciesTree, "__", 1);

        minDL = GeneSpeciesTreeUtil::Instance()->GetDLScore(tmp, speciesTree, lca_mapping_tmp);

        minNode = NULL;  //meaning, subdivide the trees
        tmp->RemoveChildren();
        delete tmp;

        //try grafting otherTree2 on every branch of otherTree
        //normally we'd use a TreeIterator but here we'll modify the tree, so better not
        vector<Node*> markedTree2Nodes = markedTree2->GetPostOrderedNodes();
        for (int imn = 0; imn < markedTree2Nodes.size(); imn++)
        {
            Node* mn = markedTree2Nodes[imn];
            if (mn != markedTree2)
            {
                Node* pmn = mn->GetParent();
                bool wasNodeDup = GeneSpeciesTreeUtil::Instance()->IsNodeDup(pmn, lca_mapping_tmp);

                Node* graftn = mn->GraftOnParentEdge(otherTree2);

                //yeah, we have to REDO this because graftn has no mapping
                unordered_map<Node*, Node*> lcamap_pnm = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(markedTree2, speciesTree, "__", 1);

                bool isNodeDup = GeneSpeciesTreeUtil::Instance()->IsNodeDup(pmn, lcamap_pnm);

                int dl = GeneSpeciesTreeUtil::Instance()->GetDLScore(markedTree2, speciesTree, lcamap_pnm);


                if (dl < minDL && !(!wasNodeDup && isNodeDup && mustRetainDupSpec))
                {

                    minDL = dl;
                    minNode = mn;
                }

                //ungraft
                graftn->RemoveChildren();
                pmn->RemoveChild(graftn);
                pmn->AddSubTree(mn);
                delete graftn;

            }
        }

        if (!minNode)
        {
            Node* tmp = new Node(false);
            tmp->AddSubTree(markedTree2);
            tmp->AddSubTree(otherTree2);
            lca_mapping[tmp] = lca_mapping[geneTree];
            return tmp;
        }
        else
        {

            minNode->GraftOnParentEdge(otherTree2);
            lca_mapping[markedTree2] = lca_mapping[geneTree];
            return markedTree2;
        }


    }
    //nothing is marked, recurse
    else
    {
        Node* t1 = minTRS_Rec(ch1, speciesTree, lca_mapping, markedGeneNodes, mustRetainDupSpec);
        Node* t2 = minTRS_Rec(ch2, speciesTree, lca_mapping, markedGeneNodes, mustRetainDupSpec);
        Node* newnode = new Node(false);
        newnode->AddSubTree(t1);
        newnode->AddSubTree(t2);

        lca_mapping[newnode] = lca_mapping[geneTree];

        return newnode;
    }

}


//geneSubtreeRoots must partition the leafset.  We don't check for it.
//NOTE : geneTree WILL BE MODIFIED AND LOSE ITS LCA MAPPING.
//Don't hold to it too closely. You have been warned.
Node* GeneSubtreeCorrector::GetSubtreeTripletRespectingHistory(Node* geneTree, Node* speciesTree,
                                         unordered_map<Node*, Node*> &lca_mapping,
                                         vector<Node *> &geneSubtreeRoots,
                                         bool mustPreserveDupSpec)
{
    if (geneSubtreeRoots.size() <= 1 || speciesTree->IsLeaf())
    {
        Node* copy = new Node(false);
        copy->CopyFrom(geneTree);
        return copy;
    }

    set<Node*> geneSubtreeRoots_set;
    for (int i = 0; i < geneSubtreeRoots.size(); i++)
    {
        geneSubtreeRoots_set.insert(geneSubtreeRoots[i]);
    }

    return minTRS_Rec(geneTree, speciesTree, lca_mapping, geneSubtreeRoots_set, mustPreserveDupSpec);


    /*if (geneSubtreeRoots.size() <= 1 || speciesTree->IsLeaf())
    {
        Node* copy = new Node(false);
        copy->CopyFrom(geneTree);
        return copy;
    }

   // vector<Node*> v;
   // v.push_back(geneTree);
   // GeneSpeciesTreeUtil::Instance()->LabelInternalNodesUniquely(v);

    //sacrifice some O(n) space for faster O(1) searching
    unordered_set<Node*> geneSubtreeRoots_set;
    for (int i = 0; i < geneSubtreeRoots.size(); i++)
    {
        geneSubtreeRoots_set.insert(geneSubtreeRoots[i]);
    }

    //call a node "marked" if it is a root of geneSubtreeRoots
    //step 1 = annotate each node with # of marked descendants
    //we WILL (previously COULD) skip this in more clever ways, but meh...


    //sometimes 2 trees get handled in one iteration, so we keep track
    //of what is handled
    unordered_set<Node*> handledSubtrees;
    vector<Node*> nodesPendingDeletion;  //we NEED to prevent deleting subtrees that are in the subtreeroots


    for (int i = 0; i < geneSubtreeRoots.size(); i++)
    {

        Node* n = geneSubtreeRoots[i];

        if (handledSubtrees.find(n) != handledSubtrees.end())
            continue;   //evil continue


        //n can't be the root (assuming input is ok)
        Node* sibl = n->GetRightSibling();
        if (!sibl)
            sibl = n->GetLeftSibling();

        Node* parent = n->GetParent();

        //since geneSubtreeRoots partition the leaves, eith sibl is marked (cherry case)
        //or it has 2+ descendants that are marked
        bool is_siblmarked = (geneSubtreeRoots_set.find(sibl) != geneSubtreeRoots_set.end());



        if (is_siblmarked)  //cherry = merge trees
        {
            handledSubtrees.insert(n);
            handledSubtrees.insert(sibl);

            Node* grandParent = NULL;
            if (!n->GetParent()->IsRoot())
                grandParent = parent->GetParent();
            Node* parentLCAMapping = lca_mapping[parent];

            vector<Node*> subtrees;
            subtrees.push_back(n);
            subtrees.push_back(sibl);
            vector< unordered_map<Node*, Node*> > subtrees_lca_mapping;



            //why twice the same line below?  because supergenetreemaker needs one lca mapping
            //per gene tree.  But in our case the same lca mapping can be used for both
            subtrees_lca_mapping.push_back(lca_mapping);
            subtrees_lca_mapping.push_back(lca_mapping);


            SuperGeneTreeMaker supermaker;

            pair<Node*, int> res = supermaker.GetSuperGeneTreeMinDL(subtrees, subtrees_lca_mapping, speciesTree, mustPreserveDupSpec);


            Node* newsubtree = res.first;

            lca_mapping[newsubtree] = parentLCAMapping;

            //now replace the parent subtree with the result.  Let's not forget to update lca mapping
            //and remember that n and sibl are handled
            if (grandParent)
            {
                grandParent->RemoveChild(parent);
                grandParent->AddSubTree(newsubtree);

                nodesPendingDeletion.push_back(parent);
            }
            else
            {

                //ah, so we end up here because there's no grandParent.  That means the geneTree root
                //is itself a cherry.  In this case, just return the result from supermaker
                return newsubtree;
            }





        }
        else    //sibling has 2+ marked descendants
        {

            if (n->IsLeaf())
                continue;   //evil continue


            Node* n1 = n->GetChild(0);
            Node* n2 = n->GetChild(1);


            //try grafting on every possible branch
            TreeIterator* it_3desc = n->GetPostOrderIterator();
            while (Node* ndesc = it_3desc->next())
            {
                if (ndesc != n)
                {
                    //graft on ndesc - parent branch and check the score

                }
            }
            n->CloseIterator(it_3desc);


            //3 options:
            //do nothing (sibl, (n1, n2))
            //move the n1 subtree on the parent-sibl branch  ((sibl, n1), n2)
            //move the n2 subtree on the parent-sibl branch  ((sibl, n2), n1)

            //option 1 costs: (sibl, (n1, n2)n)parent
            bool opt1_n_dup = (lca_mapping[n] == lca_mapping[n1] || lca_mapping[n] == lca_mapping[n2]);
            bool opt1_parent_dup = (lca_mapping[parent] == lca_mapping[n] || lca_mapping[parent] == lca_mapping[sibl]);

            //four branches can have losses: sibl-parent, n1-n, n2-n, n-parent
            int opt1_nbLosses =
                GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(lca_mapping[sibl], lca_mapping[parent], opt1_parent_dup) +
                GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(lca_mapping[n1], lca_mapping[n], opt1_n_dup) +
                GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(lca_mapping[n2], lca_mapping[n], opt1_n_dup) +
                GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(lca_mapping[n], lca_mapping[parent], opt1_parent_dup);

            int opt1_cost = opt1_nbLosses + (opt1_n_dup ? 1 : 0) + (opt1_parent_dup ? 1 : 0);


            //option 2 costs: ((sibl, n1)x, n2)xparent
            //we don't actually construct option 2...we just compute the cost it'd induce
            //the lower internal node that WOULD get created is (virtually) called x
            Node* lcamap_x = lca_mapping[sibl]->FindLCAWith(lca_mapping[n1]);
            Node* lcamap_xparent = lcamap_x->FindLCAWith(lca_mapping[n2]);
            bool opt2_x_dup = (lcamap_x == lca_mapping[sibl] || lcamap_x == lca_mapping[n1]);
            bool opt2_xparent_dup = (lcamap_xparent == lcamap_x || lcamap_xparent == lca_mapping[n2]);  //actually equal to lca_mapping[parent], but whuteva
            //four branches can have losses: sibl-x, n1-x, x-xparent, n2-xparent
            int opt2_nbLosses =
                GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(lca_mapping[sibl], lcamap_x, opt2_x_dup) +
                GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(lca_mapping[n1], lcamap_x, opt2_x_dup) +
                GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(lcamap_x, lcamap_xparent, opt2_xparent_dup) +
                GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(lca_mapping[n2], lcamap_xparent, opt2_xparent_dup);

            int opt2_cost = opt2_nbLosses + (opt2_x_dup ? 1 : 0) + (opt2_xparent_dup ? 1 : 0);

            //check that we're not messing up the n root dup/spec (if necessary)
            if (mustPreserveDupSpec)
            {
                bool was_dup = GeneSpeciesTreeUtil::Instance()->IsNodeDup(n, lca_mapping);
                if (was_dup != opt2_xparent_dup)
                    opt2_cost = 9999999;    //we don't want it
            }


            //option 3 costs: ((sibl, n2)y, n1)yparent
            //we don't actually construct option 3...we just compute the cost it'd induce
            //the lower internal node that WOULD get created is (virtually) called y
            //more or less copy pasted from above
            Node* lcamap_y = lca_mapping[sibl]->FindLCAWith(lca_mapping[n2]);
            Node* lcamap_yparent = lcamap_y->FindLCAWith(lca_mapping[n1]);
            bool opt3_y_dup = (lcamap_y == lca_mapping[sibl] || lcamap_y == lca_mapping[n2]);
            bool opt3_yparent_dup = (lcamap_yparent == lcamap_y || lcamap_yparent == lca_mapping[n1]);
            //four branches can have losses: sibl-y, n2-y, y-yparent, n1-yparent
            int opt3_nbLosses =
                GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(lca_mapping[sibl], lcamap_y, opt3_y_dup) +
                GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(lca_mapping[n2], lcamap_y, opt3_y_dup) +
                GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(lcamap_y, lcamap_yparent, opt3_yparent_dup) +
                GeneSpeciesTreeUtil::Instance()->GetNbLossesOnBranch(lca_mapping[n1], lcamap_yparent, opt3_yparent_dup);

            int opt3_cost = opt3_nbLosses + (opt3_y_dup ? 1 : 0) + (opt3_yparent_dup ? 1 : 0);

            //check that we're not messing up the n root dup/spec (if necessary)
            if (mustPreserveDupSpec)
            {
                bool was_dup = GeneSpeciesTreeUtil::Instance()->IsNodeDup(n, lca_mapping);
                if (was_dup != opt3_yparent_dup)
                    opt3_cost = 9999999;    //we don't want it
            }


            //get best scenario and apply it, if necessary
            if (opt2_cost < opt1_cost && opt2_cost <= opt3_cost)
            {
                n->RemoveChild(n1);
                sibl->GraftOnParentEdge(n1);
                //we COULD update lca mapping but we're done here, so no use.
                //also, we SHOULD delete n as it is now of degree 2...
                //but we can't, since we're iterating over n
                //So we'll delete degree 2 nodes in the end before returning
            }
            else if (opt3_cost < opt1_cost && opt3_cost < opt2_cost)
            {
                n->RemoveChild(n2);
                sibl->GraftOnParentEdge(n2);
            }

        }


    }


    for (int i = 0; i < nodesPendingDeletion.size(); i++)
    {
        delete nodesPendingDeletion[i];
    }
    nodesPendingDeletion.clear();

    Node*  gcopy = new Node(false);
    gcopy->CopyFrom(geneTree);

    gcopy->DeleteSingleChildDescendants();

    return gcopy;
    */

}
