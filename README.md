## LabelGTC
This is the implementation of the LabelGTC algorithm, a general framework for genetree correction under a duplication and loss model. 

The algorithm has been described here:

    El-Mabrouk, N., & Ouangraoua, A. (2017). A General Framework for Gene Tree Correction Based on Duplication-Loss Reconciliation. (M. Herbstritt, Ed.). Schloss Dagstuhl - Leibniz-Zentrum fuer Informatik GmbH, Wadern/Saarbruecken, Germany. https://doi.org/10.4230/lipics.wabi.2017.8



## Installation

Clone or download the github repository and run ```python setup.py build``` then ```python setup.py install```

## Running

A binary `labelgtc` is provided. Help can be obtained with ```labelgtc --help``` :
    
    usage: labelgtc [-h] -s SPECIETREE [-S SMAP] -g GENETREE [-o OUTFILE]
                    [--sep GENE_SEP] -c COVSET [--spos {prefix,postfix}] --seuil
                    SEUIL [--cost COSTDL COSTDL] [--debug]
    
    LabelGTC v1.0.1rc
    
    optional arguments:
      -h, --help            show this help message and exit
      -s SPECIETREE, --sptree SPECIETREE
                            Either the filename or the newick string of the
                            species tree.
      -S SMAP, --sMap SMAP  Gene to species map. Use the standard format.
      -g GENETREE, --gtree GENETREE
                            Either the filename or the newickg of the genetree
      -o OUTFILE, --output OUTFILE
                            Name of your output files with the corrected tree.
                            When batch is specified, each corrected genetree will
                            be printed in the appropriate output file. The
                            genetree is printed on stdout if omitted.
      --sep GENE_SEP        Gene-Specie separator for each leaf name in the
                            genetree.
      -c COVSET, --covset COVSET
                            Covering set of trees: either a list of trees
                            separated by ';' or a filename
      --spos {prefix,postfix}
                            The position of the specie name according to the
                            separator. Supported option are prefix and postfix
      --seuil SEUIL         Branch contraction threshold, when the tree is binary.
                            Use only when the tree is binary.
      --cost COSTDL COSTDL  Not implemented yet | D L : 2 float values,
                            duplication and loss cost in this order
      --debug               Debug mode
    

### Example with some test data 
```
python labelgtc -s "((A,B),C);" -g "((((a_A,x_B)0.2,(b_B,e_C)0.2)0.2,y_C)0.2,((i_B,k_A)0.1,((c_C, j_A)0.1,(d_B,(g_C,h_A)0.2)0.8)0.2)0.8)0.2;" -c "a_A;x_B;y_C;(b_B,e_C);(c_C,j_A);(g_C,h_A);d_B;(i_B,k_A);" --seuil 0.7 --debug
```
