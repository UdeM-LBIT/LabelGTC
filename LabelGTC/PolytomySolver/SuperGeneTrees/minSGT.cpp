#include "minSGT.h"



string DoSubtreeCorrection(string gcontent, string scontent, bool preserveDupSpec, string markedNodesMode, string outputmode)
{

    Node* speciesTree = NewickLex::ParseNewickString(scontent, true);
    Node* geneTree = NewickLex::ParseNewickString(gcontent, false);


    if (outputmode == "stats")
        cout<<geneTree->GetLeafSet().size()<<",";
    unordered_map<Node*, Node*> lcamap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, "__", 1);

    GeneSpeciesTreeUtil::Instance()->GetDLScore(geneTree, speciesTree, lcamap);

    int DUPS_before = GeneSpeciesTreeUtil::Instance()->LASTNBDUPS;
    int LOSSES_before = GeneSpeciesTreeUtil::Instance()->LASTNBLOSSES;

    //find subtree roots.  Their label is 'm', unless markedNodesMode = 'highspecs'
    vector<Node*> markedGeneTreeNodes;

    if (markedNodesMode == "")
    {
        TreeIterator* it = geneTree->GetPostOrderIterator();
        while (Node* n = it->next())
        {
            string lbl = n->GetLabel();
            if (lbl.length() > 0 && Util::Trim(lbl).substr(lbl.length() - 1, 1) == "m")
            {
                markedGeneTreeNodes.push_back(n);

                //remove extra m, so we don't mess gene species mapping
                if (n->IsLeaf())
                    n->SetLabel(lbl.substr(0, lbl.length() - 1));
            }
        }
        geneTree->CloseIterator(it);
    }
    else if (markedNodesMode == "highspecs")
    {
        markedGeneTreeNodes = GeneSpeciesTreeUtil::Instance()->GetGeneTreeHighestSpeciations(geneTree, speciesTree, lcamap);
    }


    if (outputmode == "stats")
    {
        cout<<markedGeneTreeNodes.size()<<",";
    }


    GeneSubtreeCorrector gsc;
    Node* newgenetree = gsc.GetSubtreeTripletRespectingHistory(geneTree, speciesTree, lcamap, markedGeneTreeNodes, preserveDupSpec);

    string result = "";

    if (!newgenetree)
    {
        cout<<"No tree was returned.  This should never happen.  Never.  The fact that you are facing the physically impossible makes us question the laws of the universe which were believed true."<<endl;
    }
    else
    {
        result = NewickLex::ToNewickString(newgenetree);

        //one last reconciliation cost
        unordered_map<Node*, Node*> lcamap_new = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(newgenetree, speciesTree, "__", 1);
        GeneSpeciesTreeUtil::Instance()->GetDLScore(newgenetree, speciesTree, lcamap_new);
        int DUPS_aftercorrection = GeneSpeciesTreeUtil::Instance()->LASTNBDUPS;
        int LOSSES_aftercorrection = GeneSpeciesTreeUtil::Instance()->LASTNBLOSSES;



        if ( outputmode == "stats" )
        {
                vector<Node*> marked_new = GeneSpeciesTreeUtil::Instance()->GetGeneTreeHighestSpeciations(newgenetree, speciesTree, lcamap_new);

                cout<<DUPS_before<<","
                    <<LOSSES_before<<","
                    <<(DUPS_before + LOSSES_before)<<","
                    <<DUPS_aftercorrection<<","
                    <<LOSSES_aftercorrection<<","
                    <<(DUPS_aftercorrection + LOSSES_aftercorrection)<<","
                    <<marked_new.size();

        }
        else
        {

            cout<<result<<endl;
        }

        //cout<<NewickLex::ToNewickString(newgenetree)<<endl;
        delete newgenetree;
    }

    delete geneTree;
    delete speciesTree;

    if (markedGeneTreeNodes.size() <= 1)
        return "";
    else
        return result;

}


string DoSuperGeneTree(string gcontent, string scontent, bool preserveDupSpec, string clades_to_preserve, string treated_trees, string outputmode)
{

    Node* speciesTree = NewickLex::ParseNewickString(scontent, true);
    /*
    cout << gcontent << "\n";
    cout << scontent << "\n";
    cout << preserveDupSpec << "\n";
    cout << outputmode << "\n";
    cout << clades_to_preserve << "\n";
    cout << treated_trees << "\n";*/

    vector<string> strGeneTrees = Util::Split(gcontent, ";", false);
    vector<Node*> geneTrees;
    vector< unordered_map<Node*, Node*> > lca_mappings;

    vector<string> strClades = Util::Split(clades_to_preserve, ";", false);
    vector<Node*> clades;

    vector<string> strTrees = Util::Split(treated_trees, ";", false);
    vector<Node*> trees;


    for (int i = 0; i < strGeneTrees.size(); i++)
    {
        string tmpnewick = Util::ReplaceAll(strGeneTrees[i], "\n", "") + ";";
        Node* geneTree = NewickLex::ParseNewickString(tmpnewick);

        unordered_map<Node*, Node*> lcamap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, "__", 1);

        geneTrees.push_back(geneTree);
        lca_mappings.push_back(lcamap);


        TreeIterator* it = geneTree->GetPostOrderIterator();
        while (Node* n = it->next())
        {
            Node* smap = lcamap[n];

            if (!smap)
            {
                cout<<n->GetLabel()<<" is not mapped!";
                return "";
            }
        }
        geneTree->CloseIterator(it);
    }

    for (int i = 0; i < strClades.size(); i++)
    {
      string tmpnewick = Util::ReplaceAll(strClades[i], "\n", "") + ";";
      Node* clade = NewickLex::ParseNewickString(tmpnewick);

      clades.push_back(clade);
    }

    for (int i = 0; i < strTrees.size(); i++)
    {
      string tmpnewick = Util::ReplaceAll(strTrees[i], "\n", "") + ";";
      Node* tree = NewickLex::ParseNewickString(tmpnewick);

      trees.push_back(tree);
    }


    SuperGeneTreeMaker sgtMaker;
    pair<Node*, int> res = sgtMaker.GetSuperGeneTreeMinDL(geneTrees, clades, trees, lca_mappings, speciesTree, preserveDupSpec, true);
    Node* superTree = res.first;
    int cost = res.second;


    string result = "";
    if (!superTree)
    {
        cout<<"It seems that no solution exists."<<endl;
    }
    else
    {
        result = NewickLex::ToNewickString(superTree);
        if (outputmode == "")
        {
            cout<<"DLCOST="<<cost<<endl;
            cout<<result<<endl;
        }
        else if (outputmode == "stats")
        {
            //stats: dup, loss, rec
            unordered_map<Node*, Node*> lcamap_new = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(superTree, speciesTree, "__", 1);
            int dl = GeneSpeciesTreeUtil::Instance()->GetDLScore(superTree, speciesTree, lcamap_new);

            cout<<GeneSpeciesTreeUtil::Instance()->LASTNBDUPS<<","
                <<GeneSpeciesTreeUtil::Instance()->LASTNBLOSSES<<","
                <<dl;
        }

        delete superTree;
    }

    for (int i = 0; i < geneTrees.size(); i++)
        delete geneTrees[i];
    delete speciesTree;

    return result;
}





void DoPolytomyCorrection(string gcontent, string scontent)
{
    Node* speciesTree = NewickLex::ParseNewickString(scontent, true);
    Node* geneTree = NewickLex::ParseNewickString(gcontent, false);


    unordered_map<Node*, Node*> lcamap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, "__", 1);


    vector<Node*> highspecs = GeneSpeciesTreeUtil::Instance()->GetGeneTreeHighestSpeciations(geneTree, speciesTree, lcamap);

    if (highspecs.size() <= 2)
    {
        //nothing to do...
        int dl = GeneSpeciesTreeUtil::Instance()->GetDLScore(geneTree, speciesTree, lcamap);

        cout<<GeneSpeciesTreeUtil::Instance()->LASTNBDUPS<<","
            <<GeneSpeciesTreeUtil::Instance()->LASTNBLOSSES<<","
            <<dl;
        delete geneTree;
        delete speciesTree;
        return;
    }

    Node* polygenetree = new Node(false);

    for (int i = 0; i < highspecs.size(); i++)
    {
        Node* ti = new Node(false);
        ti->CopyFrom(highspecs[i]);
        polygenetree->AddSubTree(ti);
    }

    unordered_map<Node*, Node*> polymap = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(polygenetree, speciesTree, "__", 1);


    Node* res = PolySolver::Instance()->SolvePolytomies(polygenetree, speciesTree, polymap);


    unordered_map<Node*, Node*> resmap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(res, speciesTree, "__", 1);
    int dl = GeneSpeciesTreeUtil::Instance()->GetDLScore(res, speciesTree, resmap);

    cout<<GeneSpeciesTreeUtil::Instance()->LASTNBDUPS<<","
        <<GeneSpeciesTreeUtil::Instance()->LASTNBLOSSES<<","
        <<dl;

    delete polygenetree;
    delete res;
    delete geneTree;
    delete speciesTree;

}


string DoSGTOnHighSpecs(string gcontent, string scontent, string clades_to_preserve, string treated_trees, bool preserveDupSpec, int maxNBLeaves, int maxNBSubtrees)
{
    Node* speciesTree = NewickLex::ParseNewickString(scontent, true);
    Node* geneTree = NewickLex::ParseNewickString(gcontent, false);


    unordered_map<Node*, Node*> lcamap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, "__", 1);


    vector<Node*> highspecs = GeneSpeciesTreeUtil::Instance()->GetGeneTreeHighestSpeciations(geneTree, speciesTree, lcamap);

    vector<string> strClades = Util::Split(clades_to_preserve, ";", false);
    vector<Node*> clades;

    vector<string> strTrees = Util::Split(treated_trees, ";", false);
    vector<Node*> trees;


    if (geneTree->GetNbLeaves() > maxNBLeaves || highspecs.size() > maxNBSubtrees)
    {
        cout<<"Skipped,Skipped,Skipped";
        delete geneTree;
        delete speciesTree;
        return gcontent;
    }

    if (highspecs.size() <= 1)
    {
        //nothing to do...
        int dl = GeneSpeciesTreeUtil::Instance()->GetDLScore(geneTree, speciesTree, lcamap);

        cout<<GeneSpeciesTreeUtil::Instance()->LASTNBDUPS<<","
            <<GeneSpeciesTreeUtil::Instance()->LASTNBLOSSES<<","
            <<dl;
        delete geneTree;
        delete speciesTree;
        return gcontent;
    }


    string gstr = "";
    //DoSGT expects a list of strings, so we'll build it for him.
    //A waste of cpu time but a gain of programmer time.
    for (int i = 0; i < highspecs.size(); i++)
    {
        gstr += NewickLex::ToNewickString(highspecs[i]);
    }

    /*for (int i = 0; i < strClades.size(); i++)
    {
      string tmpnewick = Util::ReplaceAll(strClades[i], "\n", "") + ";";
      Node* clade = NewickLex::ParseNewickString(tmpnewick);

      clades.push_back(clade);
    }

    for (int i = 0; i < strTrees.size(); i++)
    {
      string tmpnewick = Util::ReplaceAll(strTrees[i], "\n", "") + ";";
      Node* tree = NewickLex::ParseNewickString(tmpnewick);

      trees.push_back(tree);
    }*/

    string res = DoSuperGeneTree(gstr, scontent, preserveDupSpec, clades_to_preserve, treated_trees, "stats");

    delete geneTree;
    delete speciesTree;

    return res;
}
