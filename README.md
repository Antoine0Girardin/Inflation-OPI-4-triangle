# Inflation-OPI-4-triangle

This project contains the code used in the paper "Violation of the Finner inequality in the four-output triangle network". 
(https://doi.org/10.48550/arXiv.2306.05922)

We use NSI inflations to bound the correlator E2.

The main file "inflation.py" contains the code to run the inflation. The level of the inflation can be changed by setting the variables n_min and n_max between 3 and 9 with n_min<=m_max. This will fix the smaller ring and the larger ring used to constrain E2.

The files "constraints.txt" and "dense.lp" contain the details about the inflation. 

The other directories contain the value of all correlators, previously computed.

The directory "saved_results" contains all the results for n_min = 3 and all possible n_max, as well as all the results for n_min=n_max, for both maximizing and minimizing E2. Each line is a different solution found by Gurobi, each column is the value of the probability in the order given in the folder "name_probabilities". Note that all probabilities are given in the results, so the three first column are the value of the three probabilities in the triangle, then the seven next probabilities are in the square and so on. In total, the number of different probabilites for each level of the inflation is 3, 7, 11, 33, 73, 237, 703.
