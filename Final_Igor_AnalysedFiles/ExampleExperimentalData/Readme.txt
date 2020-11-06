An example igor file is provided that can be used to go through the analysis procedure. The corresponding data file is 
also provided "ExampleData.txt". The columns (called 'waves' in Igor) in the file are:

xpos_PSD : the position of the bead at fixed trap position to measure the corner frequency and Diffusion.
xpos : The position of the bead.
TrapPos: The position of the trap.
State: The state variable is a column of 0 and 1, where 0 means the ratchet protocol is on and 1 means the bead is being pulled
	back to the starting position of the trajectory.

Note that the first column (xpos_PSD) is recorded first to calibrate the trap and then the experiment for ratchet protocol is 
done that records the columns: xpos, TrapPos and State. The data is provided together for simplicity of analysis in the Igor Pro.
The position data in the first three column are in arbitrary units, multiply them by 84.8 to convert them into nanometers. The 
sampling time was 20 us as described in the main text.


==========Analysis Procedure====================

The analysis file is pre-loaded in the "Procedure" file in "Analysis_eg.pxp".

--------Trap Stiffness Analysis-----------------

1) Open the file "Analysis_eg.pxp" in Igor Pro. The four data columns and analysis file are pre-loaded.
2) In the command window run the code in Igor: PSD("xpos_PSD")
3) This produces the power spectrum density of the position of the bead.
4) Go to "Analysis" in the tool bar and click "curve fitting" that opens up the "Curve Fitting" window.
5) In the "Function and Data" window choose "Discrete_lorrentzian" from the dropdown option in "Function" box.
6) In the "Data" box choose "PS_xpos_PSD" from the dropdown.
7) Go to the next window "Data Options" and under "Weighting" choose "PSsd_xpos_PSD" from the dropdown options.
8) Go to the next window "Coefficients" and in the initial guess for the parameters input the values around:
	w_0 = 40 (corner frequency)
	w_1 = 22 (diffusion in arbitrary units)
	w_2 = 2e-9 (noise floor in arbitrary units)
9) Go to the window "Output Options" and in the box in front of "Length" change it from "Auto" to "2000". (This is optional).
10) Click "Do it". This fits the power spectrum to a Discrete Lorrentzian.
11) These produce the Fit parameters:
	w_0 = 46.77 
	w_1 = 22.66
	w_2 = 4.32 e-9

--------Ratchet output analysis-----------------

12) In the command window run the code: DoAll(xpos, trapPos, state, 46.77, 23.66)
    *The numbers 46.77 and 23.66 must be replaced by the corner frequency and diffusion values obtained in the first part of the 
    analysis.
13) This prints in the command window the "delta_g" value, the mean extrated power "P" and the input trap work "P_trap".
14) In the "Data Browser", the analysed values can be found in the waves:
	F_power : 	the final mean extracted power "P" in scaled units
	F_power_se: 	standard error of the mean of extracted power
	F_velocity:	mean velocity in units of (um/s)
	F_velocity_se:	standard error of the mean velocity
	P_trap:		mean trap work "P_trap" in scaled units
	Ptrap_se	standard error of the mean trap work
