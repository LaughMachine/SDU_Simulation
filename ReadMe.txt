The primary files associated with this readme are:
- ICUsimfunc_Clean_cmdline.py
- ICUsimfunc_Clean_B_C.py
- Sensitivity_Analysis.py
- SDUequations.py

Extras:
- ICUsimfunc_Clean_B_E.py (Used in Sensitivity_Analysis.py)

-------ICUsimfunc_Clean_cmdline.py-------
This script is run by calling:

"ICUsimfunc_Clean_cmdline.py <Sim Time> <Warm-up Time> <p_c_sc> <p_c_c_d> <p_c_sc_d> <p_sc_c> <p_sc_sc> <mu_C> <sig_C> <mu_SC> <sig_SC> <multiplier> <lb_C> <lb_S> <theta> <theta_ret> <ICU allocation> <queue length>"

Where the values in <> are the arguments taken by the script 

Example:

ICUsimfunc_Clean_cmdline.py 1001000 1000 .65 .07 .07 .07 .07 2.5 3.75 1.2 1.8 1.5 9 9 1 1 17 2


The variables used for the argument and a brief description:
Sim Time - Length of time overall simulation is run
Warm-up Time - Time simulation runs without collecting data
p_c_sc - probability patient goes from critical to semi-critical
p_c_c_d - probability patient, given that he doesn't go to semi-critical from critical, returns to critical after a delay
p_c_c_d - probability patient, given that he doesn't go to semi-critical from critical, returns to semi-critical after a delay
p_sc_c - probability patient goes from semi-critical to critical
P_sc_sc - probability patient goes from semi-critical back to semi-critical
mu_C - average serivce time in the ICU
sig_C - standard deviation of service time in ICU
mu_SC - average service time in the SDU
sig_SC - standard deviation of service time in SDU
multiplier - multiplied to service time of critical patient when being treated in SDU instead of ICU
lb_C - arrival rate of critical patients
lb_S - arrival rate of semi-critical patients
theta - average abandonment time in ICU queue
theta_ret - average time until patient returns to system, given that the patient returns
ICU allocation - number of nurses out of 20 that are allocated to the ICU
queue length - upper limit of values to simulate for the length of ICU queue line. Example: 3 means the simulation simulates the system with values of infinity and 0 to 3.


For specific of how the simulation runs, check the comments in the code, they should be comprehensive.

-------ICUsimfunc_Clean_B_C.py-------
This performs the same simulation as ICUsimfunc_Clean_cmdline.py, but the variables are written in the main function.


-------Sensitivity_Analysis.py-------
Requires ICUsimfunc_Clean_B_E.py. Takes 1 arguments: mode
Usage: "python ICUsimfunc_Clean_B_E.py <mode>"

The mode is which parameter to perform sensitivity analysis on. The code splits the values simulated into two groups in order to reduce the number of simulations run when the script is run. The mode and parameters tested are as follows:

1,2: arrival rate of semi-critical patients
3,4: rate which service time is extended when critical patients are treated in the SDU
5,6: probability return to critical
7,8: probability return to semi-critical
9,10: SDU service time
11,12: ICU service time
13,14: probability patient goes to semi-critical after critical stage
15,16: theta


-------SDUequations.py-------
Contains functions solve for optimal allocation based on equations found in the paper. The script contains commented out code demonstrating the usage of these functions. Usually, this code is to be called from other scripts. View code comments for more detail.

