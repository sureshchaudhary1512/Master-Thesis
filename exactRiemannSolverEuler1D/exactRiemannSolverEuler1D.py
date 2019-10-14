#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt
import pdb

"""
Exact (semi-)analytical solution of the Riemann problem defined with the Euler set of equations in 1D.
Based on its solution procedure explained in chap4. of Toro's book.
"""

#Define 1D geometry and mesh
nbCells = 100
#Intercells points
length = 1.0
x = np.linspace(0.,length,nbCells+1)
dx = x[1]-x[0]
y=x+(dx/2.0)
#Centres of cells
y=y[:(len(y)-1)]

#Define parameters
timeOut = 0.035
gamma = 1.4  #Ratio of specific heats
initdisc = nbCells/2  #Initial discontinuity position
posInitDisc = 0.5
#initdisc*length/dx

##Initialization of fields (Initial Conditions)
"""
#Test 1: Sod test problem. Left rarefaction, contact and right shock.
rhoL = 1.0   #Initial density on left state
vL = 0.0     #Initial velocity on left state
pL = 1.0     #Initial pressure on left state
rhoR = 0.125   #Initial density on right state
vR = 0.0     #Initial velocity on right state
pR = 0.1     #Initial pressure on right state

#Test 2: 123 problem. Two strong rarefaction waves, stationary contact discontinuity.
rhoL = 1.0   #Initial density on left state
vL = -2.0     #Initial velocity on left state
pL = 0.4     #Initial pressure on left state
rhoR = 1.0   #Initial density on right state
vR = 2.0     #Initial velocity on right state
pR = 0.4     #Initial pressure on right state

#Test 3: left half of blast problem of Woodward and Collela. 
#Left rarefaction, contact and right shock.
rhoL = 1.0   #Initial density on left state
vL = 0.0     #Initial velocity on left state
pL = 1000.0     #Initial pressure on left state
rhoR = 1.0   #Initial density on right state
vR = 0.0     #Initial velocity on right state
pR = 0.01     #Initial pressure on right state

#Test 4: right half of Woodward and Collela problem.
#left shock, contact and right rarefaction.
rhoL = 1.0   #Initial density on left state
vL = 0.0     #Initial velocity on left state
pL = 0.01     #Initial pressure on left state
rhoR = 1.0   #Initial density on right state
vR = 0.0     #Initial velocity on right state
pR = 100.0     #Initial pressure on right state
"""
#Test 5: right and left shocks
rhoL = 5.99924   #Initial density on left state
vL = 19.5975     #Initial velocity on left state
pL = 460.894     #Initial pressure on left state
rhoR = 5.99242   #Initial density on right state
vR = -6.19633     #Initial velocity on right state
pR = 46.0950     #Initial pressure on right state



rho = np.zeros(nbCells)
v = np.zeros(nbCells)
p = np.zeros(nbCells)

#Compute gamma related constants
G1 = (gamma-1.0)/(2.0*gamma)
G2 = (gamma+1.0)/(2.0*gamma)
G3 = (2.0*gamma)/(gamma-1.0) 
G4 = 2.0/(gamma-1.0)
G5 = 2.0/(gamma+1.0)
G6 = (gamma-1.0)/(gamma+1.0)
G7 = (gamma-1.0)/2.0
G8 =  gamma-1.0

#Compute sound speeds (equation (1.35))
cL = np.sqrt(gamma*pL/rhoL)
cR = np.sqrt(gamma*pR/rhoR)

#The pressure positivity is tested for (equation 4.40)
#if (G4*(cL+cR)<=(vR-vL)):
#    print "The initial data are such that vaccum is generated, program stopped"

#Definition of functions
def Compute_pv_starRegion():
    #Compute the pressure and the velocity in the star region
    p0 = guess_p()  #Give an initial guess for the iterative solution procedure of the pressure p
    dv = vR-vL
    nbIter = 20
    TOL_P = 1.0e-6
    i=0
    CHA = 1.0e3
    while ((CHA>TOL_P) and (i<nbIter)):
        FL,FLD=compute_Fs(p0,rhoL,pL,cL)         #computation of function FL (and its slope) in the exact Riemann solver (equation (4.5))
        FR,FRD=compute_Fs(p0,rhoR,pR,cR)         #computation of function FR in the exact Riemann solver (equation (4.5))
        p_new = p0 -(FL+FR+dv)/(FLD+FRD)         #Newton update with analytical slope  (equation (4.44))
        CHA = 2.0*np.abs((p_new-p0)/(p_new+p0))  #Computation of the stop criterion (equation (4.45))
        p0 = p_new
        if (p0<0.0):
            p0 = TOL_P
        i+=1
    if (i==nbIter) or (CHA>TOL_P):
        print "The pressure didnt converge"
    vS = 0.5*(vL+vR+FR-FL)                       #Compute the velocity in the star region (equation (4.9))
    print 'Number of pressure iterations: ',i
    return p_new,vS

def guess_p():
    """
    Provide a guess value for pressure pM in the star region.
    The choice is made according to adaptive Riemann solver using the PVRS, TRRS and TSRS approximate Riemann solvers.
    """
    Q_user = 2.0
    #Compute guess pressure from PVRS Riemann solver
    CUP = 0.25*(rhoL+rhoR)*(cL+cR)
    PPV = 0.5*(pL+pR)+0.5*(vL-vR)*CUP
    PPV = np.max((0.0,PPV))
    pmin = np.min((pL,pR))
    pmax = np.max((pL,pR))
    Qmax = pmax/pmin
    if ((Qmax<=Q_user) and (pmin<=PPV) and (PPV<=pmax)):
        #Select PVRS Riemann solver
        pM = PPV
    else:
        if (PPV<pmin):
            #Select Two-Rarefaction Riemann solver
            PQ = (pL/pR)**G1
            vM = (PQ*vL/cL + vR/cR + G4*(PQ-1.0))/(PQ/cL + 1.0/cR)
            PTL = 1.0 + G7*(vL-vM)/cL
            PTR = 1.0 + G7*(vM-vR)/cR
            pM = 0.5*(pL*PTL**G3 + pR*PTR**G3)
        else:
            #Select Two-Shock Riemann solver with PVRS as estimate
            GEL = np.sqrt((G5/rhoL)/(G6*pL+PPV))
            GER = np.sqrt((G5/rhoR)/(G6*pR+PPV))
            pM = (GEL*pL+GER*pR-(vR-vL))/(GEL+GER)
    return pM

def compute_Fs(p0,rhoK,pK,cK):
    #Evaluate the pressure functions FL and FR and their derivative in the exact Riemann solver
    if (p0<=pK):
        #rarefaction wave
        PRAT = p0/pK
        Fsol = G4*cK*(PRAT**G1 - 1.0)            #equation (4.6)
        FD = (1.0/(rhoK*cK))*PRAT**(-G2)         #equation (4.37) (slope)
    else:
        #Shock wave
        AK = G5/rhoK                             #equation (4.8)
        BK = G6*pK                               #equation (4.8)
        QRT = np.sqrt(AK/(BK+p0))
        Fsol = (p0 - pK)*QRT                     #equation (4.6)
        FD = (1.0 - 0.5*(p0-pK)/(BK+p0))*QRT     #equation (4.37) (slope)
    return Fsol,FD

def sample(pM,vM,S):
    """
    Sample the solution throughout the wave pattern.
    Pressure pM and velocity vM are known in the star region
    Sampling is performed in terms of the 'speed' S = X/t.
    Sample values are rho, v, p
    """
    if (S<=vM):
        #Sampling point lies to the left of the contact discontinuity
        if (pM<=pL):
            #Left rarefaction
            SHL = vL-cL        #Speed of the head of the left rarefaction wave (equation (4.55))
            if (S<=SHL):
                #Sample point is left data state
                rhoSol = rhoL
                vSol = vL
                pSol = pL
            else:
                CML = cL*(pM/pL)**G1
                STL = vM-CML   #Speed of the tail of the left rarefaction wave (equation (4.55))
                if (S>STL):
                    #Sample point is Star left state
                    rhoSol = rhoL*(pM/pL)**(1.0/gamma)           #equation (4.23)
                    vSol = vM
                    pSol = pM
                else:
                    #Sample point is inside the left fan
                    vSol = G5*(cL+G7*vL+S)                       #equation (4.56-2)
                    C = G5*(cL+G7*(vL-S))
                    rhoSol = rhoL*(C/cL)**G4                       #equation (4.56-1)
                    pSol = pL*(C/cL)**G3                         #equation (4.56-3)
        else:
            #Left shock
            PML = pM/pL
            SL = vL - cL*np.sqrt(G2*PML+G1)                      #equation (4.52) 
            if (S<=SL):
                #Sample point is left data state
                rhoSol = rhoL
                vSol = vL
                pSol = pL
            else:
                #Sample point is star left state
                rhoSol = rhoL*(PML+G6)/(PML*G6+1.0)              #equation (4.50)
                vSol = vM
                pSol = pM
    else:
        #Sampling point lies to the right of the contact discontinuity
        if (pM>pR):
            #Right shock
            PMR = pM/pR
            SR = vR + cR*np.sqrt(G2*PMR+G1)                      #equation (4.59)
            if (S>=SR):
                #Sample point is right data state
                rhoSol = rhoR
                vSol = vR
                pSol = pR
            else:
                #Sample point is star right state
                rhoSol = rhoR*(PMR+G6)/(PMR*G6+1.0)              #equation (4.57)
                vSol = vM
                pSol = pM
        else:
            #Right rarefaction
            SHR = vR+cR                                          #Speed of the head of the right rarefaction wave, equation (4.62)
            if (S>=SHR):
                #Sampled point is right data state
                rhoSol = rhoR
                vSol = vR
                pSol = pR
            else:
                CMR = cR*(pM/pR)**G1
                STR = vM+CMR                                     #Speed of the tail of the right rarefaction wave, equation (4.62)
                if (S<=STR):
                    #Sampled point is star right state
                    rhoSol = rhoR*(pM/pR)**(1.0/gamma)           #equation (4.60)
                    vSol = vM
                    pSol = pM
                else:
                    #Sampled point is inside right fan
                    vSol = G5*(-cR+G7*vR+S)                      #equation (4.63-2)
                    C = G5*(cR-G7*(vR-S))
                    rhoSol = rhoR*(C/cR)**G4                     #equation (4.63-1)
                    pSol = pR*(C/cR)**G3                         #equation (4.63-3)
    return rhoSol,vSol,pSol

####################### Start of the computation part ###########################################
#Exact solution for pressure and velocity in the star region.
pM,vM = Compute_pv_starRegion()
print 'p_star, v_star =',pM,vM
            
#Complete solution at time timeOut is found 
for i in range(nbCells):
    # Compute the 'speed' of the point y[i] at time timeOut
    #pdb.set_trace()
    S = (y[i]-posInitDisc)/timeOut
    #Solution at point (X,t) is found
    rho[i],v[i],p[i] = sample(pM,vM,S)

#Internal energy
E = p/(rho*(gamma-1.0))

#### plot ###########################################################
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')
ax1.plot(y, rho,'b-o')
ax1.grid()
ax1.set_xlabel('Position')
ax1.set_ylabel('Density')
ax2.plot(y, v,'r-o')
ax2.grid()
ax2.set_xlabel('Position')
ax2.set_ylabel('Velocity')
ax3.plot(y,p,'k-o')
ax3.grid()
ax3.set_xlabel('Position')
ax3.set_ylabel('Pressure')
ax4.plot(y,E,'g-o')
ax4.grid()
ax4.set_xlabel('Position')
ax4.set_ylabel('Internal energy')
plt.show()

