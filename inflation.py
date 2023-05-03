# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 12:40:30 2022

@author: Girardin
"""

import numpy as np
import os
import sys
import gurobipy as gp
from gurobipy import GRB
from itertools import permutations

# remove the file containing the previous constraints
try:
    os.remove('constraints.txt')
except:
    pass

#some funtion to read the data saved in the directories corr_txt, correlators, normalisation
def readinf(n):
    f_p= open('correlators/correlators{}.dat'.format(n), 'r')
    corr=np.loadtxt(f_p, dtype=int)
    f_p.close()
    return corr.tolist()

def readinfnorm(n):
    f_p= open('normalisation/norm{}.dat'.format(n), 'r')
    corr=np.loadtxt(f_p, dtype=int)
    f_p.close()
    return corr.tolist()

def readinfstr(file):
    f_p= open(file, 'r')
    corr=np.loadtxt(f_p, dtype=str)#, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})
    f_p.close()
    return corr.tolist()

def write_constr(towrite):  
    f= open('constraints.txt', 'a')
    f.write(towrite+'\n')
    f.close()
    
#some function use in the inflation
def perm_indices(p):
    '''p a string containing the indices jkl0
    return the list of all equivalent string with the jkl symbol permuted and the 0 at the same place'''
    n_outs=3
    order = np.arange(n_outs)
    indices = np.array(['j', 'k', 'l', '0'])
    ind={i:indices[i] for i in range(n_outs+1)}
    perms=[]
    for perm in permutations(order, n_outs): #all permutation of 0,1,2,3,4,5
        new_ind = {ind[perm[i]]:i for i in range(n_outs)}
        new_ind['0']=3
        # print(new_ind)
        q=''
        for j in range(len(p)):
            q+=ind[new_ind[p[j]]]
        # print(q)
        perms.append(q)
    return perms

def compute_pindex(n_poly):
    pindex = [0]*n_poly #index that gives the position of the probs for cors[n] -> vars[i+pindex[n]] with i in range(np[n])
    for i in range(len(pindex)-1):
        ind = 0
        for j in range(i+1):
            ind += nbp[j]
        pindex[i+1]=ind
    return pindex


#loading name, value of the correlators and normalisation
correlator3=readinfstr('corr_txt/corr3.txt')  
correlator4=readinfstr('corr_txt/corr4.txt') 
correlator5=readinfstr('corr_txt/corr5.txt') 
correlator6=readinfstr('corr_txt/corr6.txt') 
correlator7=readinfstr('corr_txt/corr7.txt') 
correlator8=readinfstr('corr_txt/corr8.txt') 
correlator9=readinfstr('corr_txt/corr9.txt')  

correlator3v = readinf(3)
correlator4v = readinf(4)
correlator5v = readinf(5)
correlator6v=readinf(6)
correlator7v=readinf(7)
correlator8v=readinf(8)
correlator9v=readinf(9)

cors = [correlator3, correlator4, correlator5, correlator6, correlator7, correlator8, correlator9]
corsv = [correlator3v, correlator4v, correlator5v, correlator6v, correlator7v, correlator8v, correlator9v]

norm3 = readinfnorm(3)
norm4 = readinfnorm(4)
norm5 = readinfnorm(5)
norm6 = readinfnorm(6)
norm7 = readinfnorm(7)
norm8 = readinfnorm(8)
norm9 = readinfnorm(9)


allnorm = [norm3, norm4, norm5, norm6, norm7, norm8, norm9]
nbp = []
for i in range(len(corsv)):
    nbp.append(len(corsv[i][0]))


def inflation(n_max, n_min, time_limit=0):
    """"solve the inflation for all polygon of size n_min to n_max
    time_limit : time before we stop the inflation (in seconds)"""

    n_poly = n_max-2
    n_max-=3
    n_min-=3

    model = gp.Model()
    if time_limit!=0:
        model.setParam('TimeLimit', time_limit)
        
    ## some gurobi parameter for precision
    model.params.FeasibilityTol = 1e-9
    # model.params.OptimalityTol = 1e-2
    # model.params.NumericFocus = 3
    model.params.NonConvex = 2  # 2 #to add quadratic constraints
    
    # Add variables to model
    vars = []

    # all p
    pindex=compute_pindex(n_poly)
    for i in range(n_poly):
        for j in range(nbp[i]):
            vars.append(model.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS))
    
    # constraints for normalisation of p
    for j in range(n_poly):
        expr= gp.LinExpr()
        for i in range(nbp[j]):
            expr += allnorm[j][i]*vars[i+pindex[j]]
        model.addLConstr(expr, GRB.EQUAL, 1)

    
    def addeqconstr(n,m,i,j):
        constr_str=cors[n][i]+'='+cors[n+m][j]
        # print(constr_str)
        write_constr(constr_str)
        constr_eq.append(constr_str)
        constr_eql.append(cors[n][i])
        constr_eqr.append(cors[n+m][j])
        left = gp.LinExpr()
        right = gp.LinExpr()
        for k in range(nbp[n]):
            left += corsv[n][i][k]*vars[k+pindex[n]]
        for k in range(nbp[n+m]):
            right += corsv[n+m][j][k]*vars[k+pindex[n+m]]
        model.addLConstr(left, GRB.EQUAL, right)
        
    constr_eq=[]
    constr_eql=[]
    constr_eqr=[]
    
    for n in range(n_poly-1):
        if n<=n_max and n>=n_min:
            for i in range(len(cors[n])):
                for m in range(n_poly-n-1):
                    for j in range(len(cors[n+m+1])):
                        for pos in range(len(cors[n][i])):
                            if cors[n][i][-1-pos] == '0' and cors[n][i][:-1-pos]+'0'+cors[n][i][-1-pos:] == cors[n+m+1][j]:
                                # print(cors[n][i], '=', cors[n+m+1][j])
                                constr_str=cors[n][i]+'='+cors[n+m+1][j]
                                if constr_str not in constr_eq:
                                    
                                    addeqconstr(n,m+1,i,j)
    
    def add0constr(n,i):
        constr_str=cors[n][i]+'=0'
        write_constr(constr_str)
        constr_0.append(constr_str)
        expr= gp.LinExpr()
        for j2 in range(nbp[n]):
            expr += corsv[n][i][j2]*vars[j2+pindex[n]]
        model.addLConstr(expr, GRB.EQUAL, 0)
        
    def addquadconstr(n,i,a,si,b,sj):
        write_constr(cors[n][i]+'='+cors[a][si]+'*'+cors[b][sj])
        exprl= gp.LinExpr()
        exprr1= gp.LinExpr()
        exprr2= gp.LinExpr()
        for j2 in range(nbp[n]):
            exprl += corsv[n][i][j2]*vars[j2+pindex[n]]
        for j2 in range(nbp[a]):
            exprr1 += corsv[a][si][j2]*vars[j2+pindex[a]]
        for j2 in range(nbp[b]):
            exprr2 += corsv[b][sj][j2]*vars[j2+pindex[b]]        
        model.addQConstr(exprl, GRB.EQUAL, exprr1*exprr2) 
    
    def eq_corr2(n_in, i_in):
        '''find a simpler correlator equal to the input correlator'''
        cin= cors[n_in][i_in]
        if cin in constr_eqr:
            ind = constr_eqr.index(cin)
            corr = constr_eql[ind]
            done=False
            for n in range(len(cors)):
                for i in range(len(cors[n])):
                    if cors[n][i]==corr:
                        m=n
                        j=i
                        done=True
                        break
                if done:
                    break
        else:
            m=n_in
            j=i_in
        return m,j #cors[m][j] is best corr to use
    
    def eq_corr(n_in, i_in):
        stop=False
        while not stop:
            n,i = eq_corr2(n_in, i_in)
            if n!=n_in or i!=i_in:
                n_in=n
                i_in=i
            else:
                stop=True
        return n, i
            
    
    constr_quadl=[]
    constr_quadr=[]
    constr_quad=[]
    constr_0=[]
    already_in_quad=[]
    # # no signaling constraint 'j0j0' correlator is 0 
    for n in range(n_poly):
        if n<=n_max and n>=n_min:
            for i in range(len(cors[n])):
                for j in range(len(cors[n][0])-1):
                    if cors[n][i][j] == '0' and cors[n][i][j+1] != '0':
                        for k in range(len(cors[n][i])-j-2):
                            if cors[n][i][j+2+k] == '0':
                                nonzero=False #check if we have a correlator like 'jj0jj0'
                                for a in range(6):
                                    for b in range(6):
                                        if a+b+6<=n+3: #check to look for the right size of sub correlators
                                            for si in range(len(cors[a])):
                                                for sj in range(len(cors[b])):
                                                    if any([(cors[n][i]==cors[a][si]+permut and cors[a][si][-1]=='0') for permut in perm_indices(cors[b][sj])]): #check for cases like jj0jj0=jj0kk0=jj0*jj0 
                                                        nonzero=True
                                                       
                                if nonzero:
                                    break
                                if not nonzero:
                                    constr_str=cors[n][i]+'=0'
                                    if constr_str not in constr_0:
                                        add0constr(n,i)
                                    break
                  
                                
    #check for quadratic constraints
    for n in range(n_poly):
        if n<=n_max and n>=n_min:
            for i in range(len(cors[n])):
                for j in range(len(cors[n][0])-1):
                    if cors[n][i][j] == '0' and cors[n][i][j+1] != '0':
                        for k in range(len(cors[n][i])-j-2):
                            if cors[n][i][j+2+k] == '0':
                                nonzero=False #check if we have a correlator like 'jj0jj0'
                                for a in range(6):
                                    for b in range(6):
                                        if a+b+6<=n+3: #check to look for the right size of sub correlators
                                            for si in range(len(cors[a])):
                                                for sj in range(len(cors[b])):
                                                    if any([(cors[n][i]==cors[a][si]+permut and cors[a][si][-1]=='0') for permut in perm_indices(cors[b][sj])]): #check for cases like jj0jj0=jj0kk0=jj0*jj0 
                                                        nonzero=True
                                                        a2, si2=eq_corr(a, si) #set up simpler corr
                                                        b2, sj2=eq_corr(b, sj)
                                                        
                                                        # if constr_l+'='+constr_r not in constr_quad:
                                                        #make smaller corr size of n_min
                                                        if n_min<=a2:
                                                            a3, si3=a2, si2
                                                        else:
                                                            for i2 in range(len(cors[n_min])):
                                                                if cors[n_min][i2]== cors[a2][si2] + '0'*(len(cors[n_min][i2])-len(cors[a2][si2])):
                                                                    a3, si3 = n_min, i2
                                                                    break
                                                        if n_min<=b2:
                                                            b3, sj3=b2, sj2
                                                        else:
                                                            for i2 in range(len(cors[n_min])):
                                                                if cors[n_min][i2]== cors[b2][sj2] + '0'*(len(cors[n_min][i2])-len(cors[b2][sj2])):
                                                                    b3, sj3=n_min, i2
                                                                    break
                                                                
                                                        constr_l=cors[n][i]
                                                        constr_r=cors[a3][si3]+'*'+cors[b3][sj3]
                                                                
                                                        if constr_l not in constr_quadl:
                                                            if (cors[a3][si3]+'=0' in constr_0) or (cors[b3][sj3]+'=0' in constr_0):
                                                                constr_str=cors[n][i]+'=0'
                                                                if constr_str not in constr_0:
                                                                    add0constr(n,i)
                                                            elif constr_r in constr_quadr: #make quad constr lin if already exist
                                                                already_in_quad.append(constr_l)
                                                                ind=constr_quadr.index(constr_r)
                                                                n2=len(constr_quadl[ind])-3
                                                                ind2=cors[n2].index(constr_quadl[ind])
                                                                if cors[n][i]+'='+cors[n2][ind2] not in constr_eq:
                                                                    addeqconstr(n,n2-n,i,ind2)
                                                            elif constr_l not in already_in_quad:
                                                                constr_quad.append(constr_l+'='+constr_r)
                                                                constr_quadl.append(constr_l)
                                                                constr_quadr.append(constr_r)
                                                                addquadconstr(n,i,a3,si3,b3,sj3)
                                                    
                                if nonzero:
                                    break
                                if not nonzero:
                                    constr_str=cors[n][i]+'=0'
                                    if constr_str not in constr_0:
                                        add0constr(n,i)
                                    break
    

    
    # set objectif
    obj= gp.LinExpr()
        
    #objectif = E2
    for i in range(len(cors[n_poly-1])):
        if cors[n_poly-1][i]=='jj0'+'0'*(n_poly-1):
            # print('objectif : ', cors[n_poly-1][i], corsv[n_poly-1][i])
            for j in range(len(corsv[n_poly-1][i])):
                obj+= -corsv[n_poly-1][i][j]*vars[j+pindex[n_poly-1]] #minus to maximize

    model.setObjective(obj)
    # Solve
    model.optimize() 
    
    # Write model to a file
    model.write('dense.lp')
    try:
        x= model.getAttr('X', vars)
    except:
        x=0
    
    nSolutions = model.SolCount
    allsols=[]
    if nSolutions>1:
        for i in range(nSolutions):
            model.params.SolutionNumber =i
            allsols.append(model.getAttr('Xn', vars))
    else:
        allsols.append(x)

    print('E2 max =', -model.objval)
    return x, allsols, constr_eq, constr_0, constr_quad


def print_p(n, x):
    '''print the output in a ring with the best probabilities in ALL PREVIOUS ring'''
    # p = openp('p'+str(x)+'p3') 
    
    triangle =  readinfstr('name_probabilities/triangle.txt')
    square =    readinfstr('name_probabilities/square.txt')
    pentagon =  readinfstr('name_probabilities/pentagon.txt')
    hexagon =   readinfstr('name_probabilities/hexagon.txt')
    heptagon =   readinfstr('name_probabilities/heptagon.txt')
    octogon =    readinfstr('name_probabilities/octogon.txt')
    enneagon =   readinfstr('name_probabilities/enneagon.txt')
    
    labels = [triangle, square, pentagon, hexagon, heptagon, octogon, enneagon]
    labels2 = labels[:n-1]
    labels3 = [item for sublist in labels2 for item in sublist] #flat list
    for i in range(len(x)):
        print(labels3[i], x[i])
            
if __name__=='__main__':
    #change here which ring to use in the inflation
    #only the ring from n_min to n_max will be considered
    n_min=3
    n_max=6
    x, allsols, constr_eq, constr_0, constr_quad=inflation(n_max, n_min)
    
    ## to print the probabilities of the solution found with the probability it correspond
    # print_p(n_max, x)