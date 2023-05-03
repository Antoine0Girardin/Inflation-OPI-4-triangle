# Inflation-OPI-4-triangle

This project is the code used in the paper "Around the Finner inequality in the symmetric four output triangle network".

We use NSI inflations to bound the correlator E2.

The main file "inflation.py" contains the code to run the inflation. The level of the inflation can be changed by setting the variables n_min and n_max between 3 and 9 with n_min<=m_max. This will fix the smaller ring and the larger ring used to constrain E2.

The files "constraints.txt" and "dense.lp" contain the detail about the inflation. 

The other directories contain the value of all correlators, previously computed.
