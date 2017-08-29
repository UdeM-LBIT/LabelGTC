# LabelGTC
Implementation of LabelGTC algorithm

To test it : python -m LabelGTC.tests.test1
(from the root folder)

Some little details still missing, but the recursion on minSGT seems to work properly

# You can now run the program using the binary labelgtc (This will be included in the setup.py soon)

Example of running using (test1) data :

```
python labelgtc -s "((A,B),C);" -g "((((a_A,x_B)0.2,(b_B,e_C)0.2)0.2,y_C)0.2,((i_B,k_A)0.1,((c_C, j_A)0.1,(d_B,(g_C,h_A)0.2)0.8)0.2)0.8)0.2;" -c "a_A;x_B;y_C;(b_B,e_C);(c_C,j_A);(g_C,h_A);d_B;(i_B,k_A);" --seuil 0.7 --debug
```

For now, copy the labelgtc file from LabelGTC/bin to the LabelGTC (directory with lib and src folder) before executing it.