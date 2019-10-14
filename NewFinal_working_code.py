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
from matplotlib import pyplot as plt
#Final code 


########################################### Parameters ##############################################################
rho = 1.0
C = 1.0
K = 1.0
alpha = 1.0
T0 = 200.0
a = -1/(2*rho*C)  # -ve sign
b = np.sqrt(rho*C*K)/T0


S_L =  0.1  
S_R =  0.5
B_L =  1.0
B_R =  1.0
h_L = (K*B_L/T0)
h_R = (K*B_R/T0)

##########################  Integral Curves and Hugoniot Locus ############################################################
# If S_L < S_R Left Shock and Right Rarefaction
def IntegralCurve1(S_star,S_R):
    IC1 = -(a/b)*((-2*b/alpha)*np.log((alpha - b*np.exp(S_star/a))/(alpha - b*np.exp(S_R/a))) - np.exp(-S_star/a) + np.exp(-S_R/a) + (2*b/(alpha*a))*(S_star - S_R))
    IC2 = -(2*a/alpha)*(np.log(alpha - b*np.exp(S_star/a))*(S_star/a) - np.log(alpha - b*np.exp(S_R/a))*(S_R/a) + polylog(2,b*np.exp(S_star/a)) - polylog(2,b*np.exp(S_R/a)))
    IC3 = (1/b)*((np.exp(-S_star/a))*(a - S_star) - np.exp(-S_R/a)*(a - S_R)) + ((S_star)**2 - (S_R)**2)/(a*alpha)
    B_star_I1 = B_R + IC1 + IC2 + IC3
    return B_star_I1

def Hugoniot_locus1(S_star,S_L):
    loc1 = B_L - (2*T0/np.sqrt(K))*(np.sqrt(S_L - S_star))*((np.sqrt(np.exp(S_L/(rho*C)) - np.exp(S_star/(rho*C)))))
    return loc1

# S_L > S_R Left Rarefaction and Right Shock
def IntegralCurve2(S_star,S_L):
    ic1 = -(a/b)*((-2*b/alpha)*np.log((alpha + b*np.exp(S_star/a))/(alpha + b*np.exp(S_L/a))) - np.exp(-S_star/a) + np.exp(-S_L/a) - (2*b/(alpha*a))*(S_star - S_L))
    ic2 = (2*a/alpha)*(np.log(alpha + b*np.exp(S_star/a))*(S_star/a) - np.log(alpha + b*np.exp(S_L/a))*(S_L/a) + polylog(2,-b*np.exp(S_star/a)) - polylog(2,-b*np.exp(S_L/a))) 
    ic3 = (1/b)*((np.exp(-S_star/a))*(a - S_star) - np.exp(-S_L/a)*(a - S_L)) - ((S_star)**2 + (S_L)**2)/(a*alpha)
    B_star_I2 = B_L + ic1 + ic2 + ic3
    return B_star_I2

def Hugoniot_locus2(S_star,S_R):
    loc2 = B_R + (2*T0/np.sqrt(K))*(np.sqrt(S_R - S_star))*((np.sqrt(np.exp(S_R/(rho*C)) - np.exp(S_star/(rho*C))))) 
    return loc2

#Parameter
    




Const1 = T0/np.sqrt(K)
Const2 = 2*b/alpha
Const3 = (2*b)/(a*alpha)
Const4 = (2*a)/alpha


#####################################################Compute F(S*) and dF(S*) ################################################

def Compute_F(S_star,S_L,S_R,B_L,B_R):
    
    if (S_L <= S_R):
        
 
        
        IC1 = (a/b)*(-Const2*np.log((alpha - b*np.exp(S_star/a))/(alpha - b*np.exp(S_R/a))) - np.exp(-S_star/a) + np.exp(-S_R/a) + (Const3)*(S_star - S_R))
        IC2 = (Const4)*(np.log(alpha - b*np.exp(S_star/a))*(S_star/a) - np.log(alpha - b*np.exp(S_R/a))*(S_R/a) + float(polylog(2,b*np.exp(S_star/a))) - float(polylog(2,b*np.exp(S_R/a))))
        IC3 = (1/b)*((np.exp(-S_star/a))*(a - S_star) - np.exp(-S_R/a)*(a - S_R))  
        IC4 = ((S_star)**2 - (S_R)**2)/(a*alpha)
        B_star_I1 = B_R - IC1 - IC2 + IC3 + IC4
#            #print(S_star)
#            S_star=S_L
        b_star_H1 = B_L - (2*Const1)*(np.sqrt(S_L - S_star))*(np.sqrt(np.exp(S_L/(rho*C)) - np.exp(S_star/(rho*C))))
        f_star = B_star_I1 - b_star_H1
        
        print('F(S*) =',f_star)
                 
    else: 
     
        #print('Left Rarefaction and Right Shock')
        ic1 = (a/b)*((Const2)*np.log((alpha + b*np.exp(S_star/a))/(alpha + b*np.exp(S_L/a))) - np.exp(-S_star/a) + np.exp(-S_L/a) - (Const3)*(S_star - S_L))
        ic2 = (Const4)*(np.log(alpha + b*np.exp(S_star/a))*(S_star/a) - np.log(alpha + b*np.exp(S_L/a))*(S_L/a) + float(polylog(2,-b*np.exp(S_star/a))) - float(polylog(2,-b*np.exp(S_L/a)))) 
        ic3 = (1/b)*((np.exp(-S_star/a))*(a - S_star) - np.exp(-S_L/a)*(a - S_L)) 
        ic4 = ((S_star)**2 - (S_L)**2)/(a*alpha)
        B_star_I2 = B_L - ic1 + ic2 + ic3 - ic4
        
        b_star_H2 = B_R + (2*Const1)*(np.sqrt(S_R - S_star))*(np.sqrt(np.exp(S_R/(rho*C)) - np.exp(S_star/(rho*C))))
        
        f_star = B_star_I2 - b_star_H2
        print('f(S*) =',f_star)
     
    return f_star



#=============================================================================================================================

def Compute_dF(S_star,S_L,S_R,B_L,B_R):
    if (S_L <= S_R):
     
        dF_1 =  - (a/b)*(-Const2*(-b*np.exp(S_star/a)/(a*(alpha - b*np.exp(S_star/a)))) + (np.exp(-S_star/a)/a) + Const3)
        dF_2 = - Const4*((-S_star*b*np.exp(S_star/a)/(a**2*(alpha - b*np.exp(S_star/a))) + np.log(alpha - b*np.exp(S_star/a))/a) +  float(polylog(1, b*np.exp(S_star/a))/a))
        dF_3 = (1/b)*(-np.exp(-S_star/a) - (-S_star + a)*np.exp(-S_star/a)/a)
        dF_4 = 2*S_star/(a*alpha)
        dF_5 = (2*Const1)*(-np.sqrt(np.exp(S_L/(C*rho)) - np.exp(S_star/(C*rho)))/(2*np.sqrt(S_L - S_star)) - np.sqrt(S_L - S_star)*np.exp(S_star/(C*rho))/(2*C*rho*np.sqrt(np.exp(S_L/(C*rho)) - np.exp(S_star/(C*rho)))))
        dF_star = dF_1 + dF_2 + dF_3 + dF_4 + dF_5
        print('dF(S*) =',dF_star)
        
    else:
  
        df_1 = -(a/b)*(Const2*(b*np.exp(S_star/a)/(a*(alpha + b*np.exp(S_star/a)))) + (np.exp(-S_star/a)/a) - Const3)
        df_2 = Const4*((np.log(alpha + b*np.exp(S_star/a))/a +  S_star*b*np.exp(S_star/a)/(a**2*(alpha + b*np.exp(S_star/a)))) + float(polylog(1, -b*np.exp(S_star/a))/a))
        df_3 = (1/b)*(-np.exp(-S_star/a) - (-S_star + a)*np.exp(-S_star/a)/a) - (2*S_star)/(a*alpha) 
        df_4 = (2*Const1)*(-np.sqrt(np.exp(S_R/(C*rho)) - np.exp(S_star/(C*rho)))/(2*np.sqrt(S_R - S_star)) - np.sqrt(S_R - S_star)*np.exp(S_star/(C*rho))/(2*C*rho*np.sqrt(np.exp(S_R/(C*rho)) - np.exp(S_star/(C*rho)))))
        dF_star = df_1 + df_2 + df_3 - df_4
        print('df(S*) =',dF_star)
        
    return dF_star


############################ While Loop ##############################################################################

S_star = 0.0000001    
tol = 1.0e-6
F_s1 = Compute_F(S_star,S_L,S_R,B_L,B_R)
dF_s1 = Compute_dF(S_star,S_L,S_R,B_L,B_R)
F_s = F_s1
dF_s = dF_s1
S_star_lis = []
Q_L_lis = []
Q_R_lis = []
i = 0
while(abs(F_s) > tol): 
    F_s = Compute_F(S_star,S_L,S_R,B_L,B_R)
    dF_s = Compute_dF(S_star,S_L,S_R,B_L,B_R)
    S_star_New = (S_star) - float(F_s)/(dF_s)
    S_star = S_star_New
    if ((S_L - S_star) < 0):
        break
    i+=1
    print(i)
    print('S_star = ',S_star)
    S_star_lis.append(S_star)
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
    #return S_star

#################################  Plot  ####################################################################




S = np.linspace(0.001,0.01,1000)
plLR = [] # rarefaction
plRS = [] #shock

for i in S:
    # Left Rarefaction and Right Shock
    IC1 = -(a/b)*((-2*b/alpha)*np.log((alpha - b*np.exp(i/a))/(alpha - b*np.exp(S_R/a))) - np.exp(-i/a) + np.exp(-S_R/a) + (2*b/(alpha*a))*(i - S_R))
    IC2 = -(2*a/alpha)*(np.log(alpha - b*np.exp(i/a))*(i/a) - np.log(alpha - b*np.exp(S_R/a))*(S_R/a) + polylog(2,b*np.exp(i/a)) - polylog(2,b*np.exp(S_R/a)))
    IC3 = (1/b)*((np.exp(-i/a))*(a - i) - np.exp(-S_R/a)*(a - S_R)) + ((i)**2 - (S_R)**2)/(a*alpha)
    B_star_I1 = B_R + IC1 + IC2 + IC3
    plLR.append([(B_star_I1)])

    loc1 = B_L - (2*T0/np.sqrt(K))*(np.sqrt(S_L - i))*((np.sqrt(np.exp(S_L/(rho*C)) - np.exp(i/(rho*C))))) 
    
#    print(loc1)
    plRS.append([((loc1))])
    
#print(plLR)
#print(plRS) 
#plt.plot(S,plR,label = 'U_L Rarefaction') 
#plt.plot(S,plS,label = 'U_R Shock')
#plt.plot(plLR) 
#plt.plot(plRS)
plt.legend(('U_L Rarefaction', 'U_R Shock'),loc='upper right',shadow=True)
plt.xlabel('S') 
plt.ylabel('B ')
#plt.text(0.0033, 0.992, r'S* , B* ', {'color': 'b', 'fontsize': 10})
#plt.annotate('local max', xy=(0.006,0.995), xytext=(0.1, 0.1),arrowprops=dict(facecolor='black', shrink=0.05),)
plt.title(' Left Rarefaction and Right Shock')
#plt.legend()
plt.show()







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