# Vapor-Pressure-Calculations
A program that can predict: 
(a) Ethanol vapor pressure vs. temperature 
(b) Ethanol vapor pressure vs. molar volume ethanol 
(c) The heat of vaporization vs. temperature

This program uses the Peng-Robinson equation of state.
Source used to assist in writing the program: Peng, D.Y. and Robinson, D.B., 1976. A new two-constant equation of state. Industrial & Engineering Chemistry Fundamentals, 15(1), pp.59-64.
Link: http://dns2.asia.edu.tw/~ysho/YSHO-English/2000%20Engineering/PDF/Ind%20Eng%20Che%20Fun15,%2059.pdf

Approach used for root solving: Newton-Raphson

der.m file contains the finite difference equation to find the value of the first dervivative using a five point stencil of a point x

afunction.m file contains equation 12 in the paper

To run the program go to Vapor_Pressure.m and click run in Matlab 
