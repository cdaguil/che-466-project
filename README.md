# che-466-project

Logistics

In this project you will work through the entire process of building and analyzing control systems,starting with building a process model,  then using that process model to synthesize a controllerarchitecture, and ending with comparing the performance of at least 2 different control systems.You may work in teams of 1-5 for this project; the amount of work required will be proportionalto the number of people in the team as described in the Required Work section of this document,and all members of the same team will receive the same project grade.  In the following primer,we describe the system that we want to control and list multiple tasks that you may do to buildand analyze controllers for this system.  Your work should culminate in a written final report, 1per team, which describes the critical points of your analysis of different control systems, as well asyour team’s recommendation on which to use, in a concise and cohesive story.  This project is veryopen ended - we are not seeking a “right” or “wrong” answer per se, but instead want to see youthink critically about how to approach this process.  The final report is due by 11:59PM on Friday,12/11.

Problem Statement

Cyclopentadiene (speciesA) spontaneously dimerizes to dicyclopentadiene (speciesB) in the re-versible reaction shown in Figure 1 under standard operating conditions (298 K, 1 bar).  However,high amounts ofAand low amounts ofBare required to ensure high yield and selectivity in adownstream  reactor  producing  a  rubber  copolymer  product.   To  better  ensure  this  occurs,  yourcolleagues in the design team have proposed introducing a high temperature CSTR that shifts theequilibrium towards higher concentration ofA.  The reactor has a large resistive heating elementwhich provides large quantities of heat to induce the temperature increase needed to shift the equi-librium.  Your team has been charged with designing the control system to ensure proper reactortemperature and outlet concentration ofAby manipulating the reactor input flowrate and heatingrate.  A schematic of the reactor system is shown in Figure 2.  Specific details on the reaction andreactor system are given in the ”Relevant Process Information” section.

Process Model Determination

All teams should do at least one of the following:

•Build a dynamic model for this process based on first principles mass and enery balances,resulting in a set of two nonlinear ODE’s which describe the change in reactor temperatureand concentration of A over time as a function of those two variables, the inlet flowrate andthe heating rate.  Then, linearize these ODE’s and take Laplace transforms to develop transferfunctions for the process.  You should get 4 transfer functions in all - 2 process outputs by2 process inputs.  You may make any assumptions when deriving the model that you deemappropriate, just be sure to state and justify them in your report.

•Use the supplementary Julia code to perform step tests on this system.  Use the results ofthese step tests to fit the data to transfer functions.  It is recommended that step changes areno more than±20% from the nominal steady state process inputs.  You should get 4 transferfunctions in all - 2 process outputs by 2 process inputs.  You may fit to any model type thatyou deem appropriate, just be sure to state and justify your choice in your report.

Input/Output Pairing

All teams should, based on the transfer function model(s) they derived, construct a relative gainarray  for  this  process  and  use  it  to  determine  an  appropriate  pairing  of  inputs  and  outputs  forcontrol.   Teams  should  use  this  pairing  as  the  basis  for  the  subsequent  control  system  design.Teams  may,  if  they  choose,  also  analyze  control  system  designs  based  on  the  non-recommendedinput/output pairing separately.

Control System Design and Analysis

All  teams  can  choose  between  any  of  the  following  controller  synthesis/tuning  strategies.   Notethat this list is not exhaustive, and if students are interested in a different approach, they should discuss with Prof.  Allman to determine if feasible and appropriate.  Supplementary Julia files willbe  provided  which  allow  you  to  run  the  appropriate  experiments  for  controller  tuning  and  testyour controllers’ set point tracking capabilities.  For set point tracking experiments, changes of nomore than±25 K in temperature and±0.1 M concentration from the nominal steady state processoutput values are recommended.

•Perform Ziegler-Nichols experiments on the real system to tune two SISO controllers for thesystem. You may choose which controller modes (P/PI/PID) to analyze.  Then, analyze their performance for set point tracking for both process outputs.

•Perform Cohen-Coon experiments on the real system to tune two SISO PI and PID controllersfor  the  system.   You  may  choose  which  controller  modes  (P/PI/PID)  to  analyze..   Then,analyze their performance for set point tracking for both process outputs.

•Derive two SISO ISE-near optimal controllers for this system.  Then,  analyze their perfor-mance for set point tracking for both process outputs.  Examine the controller performancefor at least two different values ofλ.

•Develop an optimization problem formulation for the model predictive control of this process.Then, analyze its performance for set point tracking for both process outputs.

Consideration of Operational Objectives

All teams should provide a brief discussion of how poor control of this process can negatively impactprocess economics, environmental and social responsibility, and safety.

Required Work

Each team ofnstudents must analyze and comparen+ 1 control systems for this process.  As an example, a solo student may choose to compare Ziegler-Nichols controllers to the ISE near-optimalcontrollers, or compare ISE near-optimal controllers designed for the linearized process model vs.those  designed  for  the  data-driven process model, or consider Ziegler-Nichols controllers for the two different input/output pairing options. A team of 5 may choose to, for example,  comparemodel predictive controllers using both the data driven and linearized models (2 systems) to ISEnear-optimal controllers for both input/output pairing options and both data-driven and linearizedmodels (4 more systems, for a total of 6).  If you are unsure if the number of different systems isappropriate for the number of people on the team, please contact Prof. Allman to ask.

The list of different options for control system design and analysis is not exhaustive - teams shoulddiscuss with Prof.  Allman if there is a different type of control system they are interested in look-ing at to determine if feasible and appropriate.  All teams must consider the input/output pairingproblem and discuss operational objectives.  All teams must submit their work in a written reportwhich  tells  a  concise  and  cohesive  story,  which  includes  at  the  end  a  recommendation  of  whichcontrol system to use out of the ones they’ve analyzed based on the analysis of the report.

Relevant Process Information

•Enthalpy of forward reaction:  -77 kJ/mol B produced

•Both forward and reverse reaction follow an Arrhenius rate law dependence, krxn=A0exp(−EA/RT)

•Forward activation energy:  63 kJ/mol•Reverse activation energy:  140 kJ/mol

•Pre-exponential term for forward rxn:  1.28*102M−1s−1

•Pre-exponential term for reverse rxn:  3.44*109s−1

•Forward reaction is second order in A, reverse reaction is first order in B.

•Density of mixture:  850 g/L, may be assumed constant•Heat capacity of process fluid:  1.74 J/g/K, may be assumed constant

•Volume of fluid in reactor:  100 L

•Molecular weight of A: 66.1 g/mol A

•Molecular weight of B: 132.2 g/mol B

•A good, fast flow control loop is already in place to manipulate the outlet flow and keep thereactor volume constant.

•Nominal inlet flow of reactor:  0.1 L/s

•Nominal heating rate of reactor heating element:  42.6 kW

•Nominal temperature of reactor:  625 K

•Nominal outlet concentration of A: 1.396 M

•Inlet concentration ofBis 1 M, inlet concentration of A is negligible.

•Inlet temperature is 373 K.

•All actuator and sensor dynamics are negligible
