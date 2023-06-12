# Inflation-OPI-4-triangle

This project contains the code used in the paper "Violation of the Finner inequality in the four-output triangle network". 
(https://doi.org/10.48550/arXiv.2306.05922)

We use NSI inflations to bound the correlator E2.

The main file "inflation.py" contains the code to run the inflation. The level of the inflation can be changed by setting the variables n_min and n_max between 3 and 9 with n_min<=m_max. This will fix the smaller ring and the larger ring used to constrain E2.

The files "constraints.txt" and "dense.lp" contain the details about the inflation. 

The other directories contain the value of all correlators, previously computed.

The directory "saved_results" contains all the results for n_min = 3 and all possible n_max, as well as all the results for n_min=n_max, for both maximizing and minimizing E2.
