

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 22:13:29 2019

@author: sureshchaudhary
"""

import numpy as np 
#from sympy import *
#import pdb
from mpmath import polylog 
#from mmpmath import *
#from scipy.lib.six import callable
#from scipy.optimize import cobyla
#from .optimize import OptimizeResult, _check_unknown_options
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
from scipy.optimize import Bounds
#from scipy import optimize
#from scipy.optimize import _cobyla
#from matplotlib import pyplot as plt
#Final code 



########################################### Parameters ##############################################################



rho = 970.0
C = 0.052
K = 875.0
alpha = 1.0
T0 = 200.0
a = -1/(2*rho*C)  # -ve sign
b = np.sqrt(rho*C*K)/T0

Const1 = T0/np.sqrt(K)
Const2 = 2*b/alpha
Const3 = (2*b)/(a*alpha)
Const4 = (2*a)/alpha
S_L =  0.1  
S_R =  0.2
B_L =  1.0
B_R =  1.0
h_L = (K*B_L/T0)
h_R = (K*B_R/T0)

#=================================== # Integral Curves and Hugoniot Locus # =====================================================================================================

#================================== #S_L < S_R Left Shock and Right Rarefaction =================================================================================================
def IntegralCurve1(S_star,S_R):
    IC1 = (a/b)*(-Const2*np.log((alpha - b*np.exp(S_star/a))/(alpha - b*np.exp(S_R/a))) - np.exp(-S_star/a) + np.exp(-S_R/a) + (Const3)*(S_star - S_R))
    IC2 = - (Const4)*(np.log(alpha - b*np.exp(S_star/a))*(S_star/a) - np.log(alpha - b*np.exp(S_R/a))*(S_R/a) + polylog(2,b*np.exp(S_star/a)) - polylog(2,b*np.exp(S_R/a)))
    IC3 = (1/b)*((np.exp(-S_star/a))*(a - S_star) - np.exp(-S_R/a)*(a - S_R)) + ((S_star)**2 - (S_R)**2)/(a*alpha)
    IC  =  IC1 + IC2 + IC3
    B_star_I1 = B_R - IC
    
    return B_star_I1
 
def Hugoniot_locus1(S_star,S_L):
     b_star_H1 = B_L - (2*Const1)*(np.sqrt(S_L - S_star))*(np.sqrt(np.exp(S_L/(rho*C)) - np.exp(S_star/(rho*C))))
     return b_star_H1

#============================================ # S_L > S_R Left Rarefaction and Right Shock=======================================================================================
def IntegralCurve2(S_star,S_L):
    ic1 = (a/b)*((Const2)*np.log((alpha + b*np.exp(S_star/a))/(alpha + b*np.exp(S_L/a))) - np.exp(-S_star/a) + np.exp(-S_L/a) - (Const3)*(S_star - S_L))
    ic2 = (Const4)*(np.log(alpha + b*np.exp(S_star/a))*(S_star/a) - np.log(alpha + b*np.exp(S_L/a))*(S_L/a) + polylog(2,-b*np.exp(S_star/a)) - polylog(2,-b*np.exp(S_L/a))) 
    ic3 = (1/b)*((np.exp(-S_star/a))*(a - S_star) - np.exp(-S_L/a)*(a - S_L)) - ((S_star)**2 - (S_L)**2)/(a*alpha)
    ic  = ic1 + ic2 + ic3
    B_star_I2 = B_L + ic
    return B_star_I2
 
def Hugoniot_locus2(S_star,S_R):
    b_star_H2 = B_R + (2*Const1)*(np.sqrt(S_R - S_star))*(np.sqrt(np.exp(S_R/(rho*C)) - np.exp(S_star/(rho*C))))
    return b_star_H2
#=================================================================================================================================================================================





#==================================================== Compute F(S*) and dF(S*) =============================================================================================================================
def Compute_F(S_star,S_L,S_R,B_L,B_R):
    
    if (S_L <= S_R): 
        
        #print('Left Shock and Right Rarefaction')
        
        IC1 = (a/b)*(-Const2*np.log((alpha - b*np.exp(S_star/a))/(alpha - b*np.exp(S_R/a))) - np.exp(-S_star/a) + np.exp(-S_R/a) + (Const3)*(S_star - S_R))
        IC2 = -(Const4)*(np.log(alpha - b*np.exp(S_star/a))*(S_star/a) - np.log(alpha - b*np.exp(S_R/a))*(S_R/a) + polylog(2,b*np.exp(S_star/a)) - polylog(2,b*np.exp(S_R/a)))
        IC3 = (1/b)*((np.exp(-S_star/a))*(a - S_star) - np.exp(-S_R/a)*(a - S_R)) + ((S_star)**2 - (S_R)**2)/(a*alpha)
        IC  =  IC1 + IC2 + IC3
        f_star = B_R - IC -  B_L + (2*Const1)*(np.sqrt(S_L - S_star))*((np.sqrt(np.exp(S_L/(rho*C)))) - (np.sqrt(np.exp(S_star/(rho*C)))))
        F_star = f_star
    
        print('F_star = ',F_star)
                 
    else: 
        
         
        #print('Left Rarefaction and Right Shock')
        ic1 = (a/b)*((Const2)*np.log((alpha + b*np.exp(S_star/a))/(alpha + b*np.exp(S_L/a))) - np.exp(-S_star/a) + np.exp(-S_L/a) - (Const3)*(S_star - S_L))
        ic2 = (Const4)*(np.log(alpha + b*np.exp(S_star/a))*(S_star/a) - np.log(alpha + b*np.exp(S_L/a))*(S_L/a) + polylog(2,-b*np.exp(S_star/a)) - polylog(2,-b*np.exp(S_L/a))) 
        ic3 = (1/b)*((np.exp(-S_star/a))*(a - S_star) - np.exp(-S_L/a)*(a - S_L)) - ((S_star)**2 - (S_L)**2)/(a*alpha)
        ic  = ic1 + ic2 + ic3
        f_star = B_L - ic - B_R - (2*Const1)*(np.sqrt(S_R - S_star))*((np.sqrt(np.exp(S_R/(rho*C)))) - (np.sqrt(np.exp(S_star/(rho*C)))))
        F_star = f_star
        
        print('F_star = ',F_star)
     
    return F_star



#=================================================================================================================================================================================

def Compute_dF(S_star,S_L,S_R,B_L,B_R):
    if (S_L <= S_R):
        
        dF_1 =  - (a/b)*(-Const2*(-b*np.exp(S_star/a)/(a*(alpha - b*np.exp(S_star/a)))) + (np.exp(-S_star/a)/a) + Const3)
        dF_2 =  - (Const4*(-S_star*b*np.exp(S_star/a)/(a**2*(alpha - b*np.exp(S_star/a))) +  polylog(1, b*np.exp(S_star/a))/a))
        dF_3 =   (1/b)*(-np.exp(-S_star/a) - (-S_star + a)*np.exp(-S_star/a)/a) + (2*S_star/(a*alpha)) 
        dF_4 =   (2*Const1)*(-np.sqrt(np.exp(S_L/(C*rho)) - np.exp(S_star/(C*rho)))/(2*np.sqrt(S_L - S_star)) - np.sqrt(S_L - S_star)*np.exp(S_star/(C*rho))/(2*C*rho*np.sqrt(np.exp(S_L/(C*rho)) - np.exp(S_star/(C*rho)))))
        dF_star = dF_1 + dF_2 + dF_3 + dF_4
        print('dF_star =',dF_star)
        
    else:
        
        df_1 = -(a/b)*(Const2*(b*np.exp(S_star/a)/(a*(alpha + b*np.exp(S_star/a)))) + (np.exp(-S_star/a)/a) - Const3)
        df_2 = Const4*(np.log(alpha + b*np.exp(S_star/a))/a +  S_star*b*np.exp(S_star/a)/(a**2*(alpha + b*np.exp(S_star/a))) + polylog(1, -b*np.exp(S_star/a))/a)
        df_3 = (1/b)*(np.exp(S_star/a)*S_star/a - 2*np.exp(-S_star/a)) + (2*S_star)/(a*alpha)
        df_4 = (2*Const1)*(-np.sqrt(np.exp(S_R/(C*rho)) - np.exp(S_star/(C*rho)))/(2*np.sqrt(S_R - S_star)) - np.sqrt(S_R - S_star)*np.exp(S_star/(C*rho))/(2*C*rho*np.sqrt(np.exp(S_L/(C*rho)) - np.exp(S_star/(C*rho)))))
        dF_star = df_1 + df_2 + df_3 + df_4
        print('dF_star =',dF_star)
        
    return dF_star



        
#===============================================  Constraint ==================================================================================================================================        
#S_L < S_R
ineq_cons1 = ({'type': 'ineq', 'fun': lambda S_star: 1 - b*np.exp(S_star/a) , 'jac': lambda S_star: -(b/a)*np.exp(S_star/a)},
              {'type': 'ineq', 'fun': lambda S_R: 1 - b*np.exp(S_R/a)},
              {'type': 'ineq', 'fun': lambda S_star: S_L-S_star},
              {'type': 'ineq', 'fun': lambda S_star: S_R-S_star},
              {'type': 'ineq', 'fun': lambda S_star: np.exp(S_L/rho*C)-np.exp(S_star/rho*C)},
              {'type': 'ineq', 'fun': lambda S_star: np.exp(S_R/rho*C)-np.exp(S_star/rho*C)})






##S_L > S_R            
#ineq_cons2 = ({'type': 'ineq', 'fun': lambda S_star: 1 + b*np.exp(S_star/a), 'jac': lambda S_star: (b/a)*np.exp(S_star/a)},
#              {'type': 'ineq', 'fun': lambda S_L: 1 + b*np.exp(S_L/a)},
#              {'type': 'ineq', 'fun': lambda S_star: S_R-S_star},
#              {'type': 'ineq', 'fun': lambda S_star: S_L-S_star},
#              {'type': 'ineq', 'fun': lambda S_star: np.exp(S_R/rho*C)-np.exp(S_star/rho*C)},
#              {'type': 'ineq', 'fun': lambda S_star: np.exp(S_L/rho*C)-np.exp(S_star/rho*C)})


################################################ Minimization ##################################################################
#S_star = 0.001
#linear_constraint = LinearConstraint([b*np.exp(S_star/a)], [0], [1])

#Sol = minimize(Compute_F,S_star,args=(S_L,S_R,B_L,B_R),method='SLSQP',jac=Compute_dF,constraints=[linear_constraint],options={'ftol': 1e-8, 'disp': True})

Sol = minimize(Compute_F,S_star,args=(S_L,S_R,B_L,B_R),method='SLSQP',jac=Compute_dF,constraints=ineq_cons1,options={'ftol': 1e-8, 'disp': True})
#print(Sol)                  
                  
        
#S_star = 0.00001  
#Sol = minimize(Compute_F,S_star,args=(S_L,S_R,B_L,B_R),method='COBYLA',constraints=ineq_cons1,tol=None, callback=None, options={'rhobeg': 1.0, 'maxiter': 1000, 'disp': False, 'catol': 0.0002})                
                  
                   
            
 










    




#def Compute_dF(S_star,S_L,S_R,B_L,B_R):
#    if (S_L <= S_R): 
#        #print('Left Shock and Right Rarefaction')
#        
#        dF_star = - (a/b)*(-Const2*(-b*np.exp(S_star/a)/(a*(alpha - b*np.exp(S_star/a)))) + (np.exp(-S_star/a)/a) + Const3) \
#                  - (Const4*(-S_star*b*np.exp(S_star/a)/(a**2*(alpha - b*np.exp(S_star/a))) +  polylog(1, b*exp(S_star/a))/a) \
#                  + (1/b)*(-np.exp(-S_star/a) - (-S_star + a)*np.exp(-S_star/a)/a) + (2*S_star/(a*alpha)) \
#                  +  (2*Const1)*(-np.sqrt(np.exp(S_L/(C*rho)) - np.exp(S_star/(C*rho)))/(2*np.sqrt(S_L - S_star)) - np.sqrt(S_L - S_star)*np.exp(S_star/(C*rho))/(2*C*rho*np.sqrt(np.exp(S_L/(C*rho)) - np.exp(S_star/(C*rho)))))
#                  
#                  
#        print('dF_star =',dF_star)
#        
#    else:          
#        #print('Left Rarefaction and Right Shock')
#        dF_star =  -(a/b)*(Const2*(b*np.exp(S_star/a)/(a*(alpha + b*np.exp(S_star/a)))) + (np.exp(-S_star/a)/a) - Const3) \
#                   + Const4*(np.log(alpha + b*np.exp(S_star/a))/a +  S_star*b*np.exp(S_star/a)/(a**2*(alpha + b*np.exp(S_star/a))) + polylog(1, -b*np.exp(S_star/a))/a) \
#                   + (1/b)*(np.exp(S_star/a)*S_star/a - 2*np.exp(-S_star/a)) + (2*S_star)/(a*alpha) \
#                   + (2*Const1)*(-np.sqrt(np.exp(S_R/(C*rho)) - np.exp(S_star/(C*rho)))/(2*np.sqrt(S_R - S_star)) - np.sqrt(S_R - S_star)*np.exp(S_star/(C*rho))/(2*C*rho*np.sqrt(np.exp(S_L/(C*rho)) - np.exp(S_star/(C*rho)))))
#                   
# 
#        
#        print('dF_star =',dF_star)
#    
#    return dF_star


###############################################Constraint###################################################################

#S_L < S_R
#ineq_cons1 = ({'type': 'ineq', 'fun': lambda S_star: 1 - b*np.exp(S_star/a), 'jac': lambda S_star: -(b/a)*np.exp(S_star/a)},
#              {'type': 'ineq', 'fun': lambda S_R: 1 - b*np.exp(S_R/a)},
#              {'type': 'ineq', 'fun': lambda S_star: S_L-S_star},
#              {'type': 'ineq', 'fun': lambda S_star: S_R-S_star},
#              {'type': 'ineq', 'fun': lambda S_star: np.exp(S_L/rho*C)-np.exp(S_star/rho*C)},
#              {'type': 'ineq', 'fun': lambda S_star: np.exp(S_R/rho*C)-np.exp(S_star/rho*C)})

##S_L > S_R            
#ineq_cons2 = ({'type': 'ineq', 'fun': lambda S_star: 1 + b*np.exp(S_star/a), 'jac': lambda S_star: (b/a)*np.exp(S_star/a)},
#              {'type': 'ineq', 'fun': lambda S_L: 1 + b*np.exp(S_L/a)},
#              {'type': 'ineq', 'fun': lambda S_star: S_R-S_star},
#              {'type': 'ineq', 'fun': lambda S_star: S_L-S_star},
#              {'type': 'ineq', 'fun': lambda S_star: np.exp(S_R/rho*C)-np.exp(S_star/rho*C)},
#              {'type': 'ineq', 'fun': lambda S_star: np.exp(S_L/rho*C)-np.exp(S_star/rho*C)})


################################################ Minimization ##################################################################
#S_star = 0.0001
#Sol = minimize(Compute_F,S_star,args=(S_L,S_R,B_L,B_R),method='SLSQP',jac=Compute_dF ,constraints=ineq_cons1,options={'ftol': 1e-8, 'disp': True})
#print(Sol)

#S_star = 0.1
#Sol = minimize(Compute_F,S_star,args=(S_L,S_R,B_L,B_R),method='COBYLA' ,constraints=ineq_cons1,options={'ftol': 1e-8, 'disp': True})










#a_star = S_star/(rho*C)
#a_L = S_L/(rho*C)
#a_R = S_R/(rho*C)
#gamma1 = np.sqrt(np.exp(a_L)- np.exp(a_star))
#gamma2 = np.sqrt(np.exp(a_R)- np.exp(a_star))
#PLog2_L = b*np.exp(S_L/a)
#PLog2_R = b*np.exp(S_R/a)
#PLog2_star = b*np.exp(S_star/a)
#dPLog2_star = -(b/a)*np.exp(S_star/a)
#Const1 = T0/np.sqrt(K)
#Const2 = 2*b/alpha
#Const3 = 2*a/(alpha*b)
#Const4 = (2*b)/(a*alpha)
#Zeta_star = alpha - PLog2_star
#Zeta_R = alpha - PLog2_R 
#Eta_star = alpha + PLog2_star
#Eta_L = alpha + PLog2_L
#
#
######################################################Compute F(S*) and dF(S*) ################################################
#
#def Compute_F(S_star,S_L,S_R,B_L,B_R):
#    if (S_L <= S_R): 
#        #print('Left Shock and Right Rarefaction')
#        F_star = B_L - 2*Const1*np.sqrt(S_L - S_star)*gamma1 - B_R + (a/b)*((-Const2)*np.log(Zeta_star/Zeta_R) - np.exp(-S_star/a) + np.exp(-S_R/a) + Const4*(S_star-S_R)) \
#                + Const3*(np.log(Zeta_star)*(S_star/a) - np.log(Zeta_R)*(S_R/a) + polylog(2,PLog2_star) - polylog(2,PLog2_R)) - (1/b)*(np.exp(-S_star/a)*(a-S_star) - np.exp(-S_R/a)*(a-S_R)) + (S_R**2 - S_star**2/(a*alpha))
#        print('F_star = ',F_star)
#                 
#    else:          
#        #print('Left Rarefaction and Right Shock')
#        F_star = B_R + 2*Const1*np.sqrt(S_R - S_star)*gamma2 - B_L + (a/b)*(Const2*np.log(Eta_star/Eta_L) - np.exp(-S_star/a) + np.exp(-S_L/a) - Const4*(S_star-S_L)) \
#                - Const3*(np.log(Eta_star)*(S_star/a) - np.log(Eta_L)*(S_L/a) + polylog(2,-PLog2_star) - polylog(2,-PLog2_L)) + (1/b)*(np.exp(-S_star/a)*(a-S_star) - np.exp(-S_L/a)*(a-S_L)) + (S_star**2 - S_L**2)/(a*alpha)
#     
#    return F_star
#
#
#
#def Compute_dF(S_star,S_L,S_R,B_L,B_R):
#                
#    if (S_L <= S_R): 
#        #print('Left Shock and Right Rarefaction')
#        dF_star = Const1*(np.sqrt(S_L-S_star)*np.exp(a_star)/(rho*C*gamma1) + gamma1/np.sqrt(S_L-S_star)) +(a/b)*((b*Const4*np.exp(S_star/a))/Zeta_star + np.exp(-S_star/a)/a + Const4) \
#                 + Const3*(np.log(Zeta_star)/a - (S_star*b*np.exp(S_star/a))/(a**2*Zeta_star) - np.log(Zeta_star)/a)  \
#                  -(1/b)*((S_star/a)*np.exp(-S_star/a) - 2*np.exp(-S_star/a)) - (2*S_star)/(a*alpha)
#        print('dF_star =',dF_star)
#    else:          
#        #print('Left Rarefaction and Right Shock')
#        dF_star = -Const1*(np.sqrt(S_R-S_star)*np.exp(a_star)/(rho*C*gamma2) + gamma2/np.sqrt(S_R-S_star)) + (a/b)*(b*Const4*np.exp(S_star/a)/Eta_star + np.exp(-S_star/a)/a - Const4) \
#                + Const3*((np.log(Eta_star)/a) + (S_star*b*np.exp(S_star/a)/((a**2)*Eta_star)) - (np.log(Eta_star)/a)) \
#                + (1/b)*((S_star*np.exp(-S_star/a))/a - 2*np.exp(-S_star/a)) + (2*S_star)/(a*alpha)  \
#    
#    return dF_star


























#S_L < S_R
#ineq_cons1 = ({'type': 'ineq', 'fun': lambda S_star: 1 - PLog2_star, 'jac': lambda S_star: dPLog2_star},
#             {'type': 'ineq', 'fun': lambda S_R: 1 - PLog2_R},
#             {'type': 'ineq', 'fun': lambda S_star: S_L-S_star},
#             {'type': 'ineq', 'fun': lambda S_star: gamma1})

#S_L > S_R            
#ineq_cons2 = ({'type': 'ineq', 'fun': lambda S_star: 1 + PLog2_star, 'jac': lambda S_star: -dPLog2_star},
#             {'type': 'ineq', 'fun': lambda S_L: 1 + PLog2_L},
#             {'type': 'ineq', 'fun': lambda S_star: S_R-S_star},
#             {'type': 'ineq', 'fun': lambda S_star: gamma2})



#sol = optimize.newton(Compute_F,0.0001,fprime=Compute_dF,args=(S_L,S_R,B_L,B_R),tol =1e-8)

############################ While Loop ##############################################################################
    
#tol = 1.0e-6
#F_s1,dF_s1 = Compute_F_dF(S_L,S_R,B_L,B_R,S_star)
#F_s = F_s1
#dF = dF_s1
#S_star_lis = []
#Q_L_lis = []
#Q_R_lis = []
#i = 0
#while(abs(F_s) > tol): 
#    F_s,dF_s = Compute_F_dF(S_L,S_R,B_L,B_R,S_star)
#    S_star_New = (S_star) - float(F_s)/(dF_s)
#    S_star = S_star_New
#    if ((S_L - S_star) < 0) and Zeta_star < 1:
#        break
#    
#    i+=1
#    print(i)
#    print('S_star = ',S_star)
#    S_star_lis.append(S_star)
#    if (S_L <= S_R): 
#        Q_L = Hugoniot_locus1(S_star,S_L) #left shock
#        Q_R = IntegralCurve1(S_star,S_R) # Right Rarefaction
#    else: 
#        Q_L = Hugoniot_locus2(S_star,S_R) # Right Shock
#        Q_R = IntegralCurve2(S_star,S_L) # left Rarefaction
#    Q_L_lis.append(Q_L)
#    print('Q_L = ',Q_L)
#    Q_R_lis.append(Q_R)
#    print('Q_R = ',Q_R)
#    #return S_star

#################################  Plot  ####################################################################

#plt.subplots_adjust(hspace=0.2)
#plt.subplot(221)
#plt.plot(Q_R_lis,S_star_lis, 'b-')
#plt.xlabel('Q_R')
#plt.ylabel('S*')
#plt.show()
#
#
#plt.subplots_adjust(hspace=0.2)
#plt.subplot(222)
#plt.plot(Q_L_lis,S_star_lis, 'r-')
#plt.xlabel('Q_L')
#plt.ylabel('S*')
#plt.show()
#
#
#plt.plot(Q_R_lis,S_star_lis, 'b-',Q_L_lis,S_star_lis, 'r-')
#plt.ylabel('S*')
#plt.xlabel('Q_R ,Q_L')
#plt.show()