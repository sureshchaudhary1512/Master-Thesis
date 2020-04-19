# Master Thesis 

# Riemann solver for non-linear Green-Naghdi II Thermal Equations

This thesis presents the Riemann solution for non linear Green Naghdi dissipationless energy thermal equation which is also know as type II. The classical theory of heat conduction based on the Fourier law allows infinite diffusion speed which is not well accepted from a physical point of view. The Green and Naghdi employ a procedure which differs from the usual one.They introduce the notion of thermal displacement so that thermal waves propagate at finite speed and constitutive equation for entropy flux is determined by potential function. This work is focused on the solution of a Riemann problem for Green-Naghdi II equations, which are hyperbolic and these equations are formulated according to some internal energy or free energy potential describing the constitutive response. The solution to the Riemann problem is derived according to a linearized thermal response as well as a nonlinear . It is shown that the solution of such a Riemann problem can be used in a finite volume solver.


# Test Case 

<img src= "Test_case1.png" >

# Results FEM vs FVM

<img src= "Test_case12.png" >

<img src= "Test_case13.png" >


# Conclusion

In this work, a Riemann solution has been proposed to non linear Green Naghdi type II thermal equations. The nonlinear thermal response extracted from the classical generalized thermoelasticity, written in the form of free energy potential, the dissipationless Green- Naghdi equations consist of a set of nonlinear hyperbolic equations. In a 1D medium, these equations consist of two genuinely nonlinear characteristic fields. It has been shown that the characteristic structure exhibited in the solution of a Riemann problem necessarily consists of a combination of a shock and a rarefaction(i.e. 1-shock 2-rarefaction or 1-rarefaction 2-shock). It has also been shown that equations linked to a shock and a rarefaction can be derived analytically, though the solution of the star state requires to find numerically the root of a nonlinear scalar equation. Finally, different test cases have been shown and integral curves plotted for a given set of material parameters associated with Bismuth.
