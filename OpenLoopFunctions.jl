#Functions for generating step test plots for open loop changes in process inputs
#Inlet flowrate and heating rate
#Last updated by AA - 19 November 2020

using Polynomials, Plots, DifferentialEquations

function loadProcessData()
    global DH_fwd=-77000 #J/mol B
    global EA_fwd=63000 #J/mol
    global EA_rev=140000 #J/mol
    global R_gas=8.31 #J/mol/K
    global A_fwd=1.28e2 #L/mol/s
    global A_rev=3.44e9 #/s
    global rho=850 #g/L
    global c_p=1.74 #J/g/K
    global V_rxtr=100 #L
    global MW_A=66.1 #g/mol A
    global MW_B=132.2 #g/mol B
    global flow_nom=0.1 #L/s
    global T_nom=625 #K
    global CB_in=1 #M
    global T_in=373 #K
    a0=2*V_rxtr*A_rev*exp(-EA_rev/R_gas/T_nom)*CB_in
    a1=-flow_nom-V_rxtr*A_rev*exp(-EA_rev/R_gas/T_nom)
    a2=-2*V_rxtr*A_fwd*exp(-EA_fwd/R_gas/T_nom)
    global CA_nom=roots(Polynomial([a0,a1,a2]))[2]
    global heat_nom=rho*c_p*flow_nom*(T_nom-T_in)+DH_fwd*V_rxtr*A_fwd*exp(-EA_fwd/R_gas/T_nom)*
        CA_nom^2-DH_fwd*V_rxtr*A_rev*exp(-EA_rev/R_gas/T_nom)*(CB_in-CA_nom/2)
    println("Parameters Loaded!")
end



function openLoopStepFlow(M_step;tmax=10000.0)
    f(y,p,t)=[-(flow_nom+M_step)*y[1]/V_rxtr-2*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        2*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2),
        (flow_nom+M_step)/V_rxtr*(T_in-y[2])-DH_fwd/rho/c_p*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        DH_fwd/rho/c_p*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2)+heat_nom/rho/V_rxtr/c_p]
    prob=ODEProblem(f,[CA_nom,T_nom],(0.0,tmax))
    soln=solve(prob,Rosenbrock23())
    t=soln.t
    CA_soln=zeros(size(soln.u)[1])
    T_soln=zeros(size(soln.u)[1])
    for i=1:size(soln.u)[1]
        CA_soln[i]=soln.u[i][1]
        T_soln[i]=soln.u[i][2]
    end
    plot(t,[CA_soln,T_soln],layout=(2,1),legend=false,xlabel="Time (s)", ylabel=["Conc. of A (M)" "Temp. (K)"],gridalpha=0.75,minorgrid=true,minorgridalpha=0.3)
#    return t,CA_soln,T_soln
end

function openLoopStepHeat(M_step;tmax=10000.0)
    f(y,p,t)=[-flow_nom*y[1]/V_rxtr-2*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        2*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2),
        flow_nom/V_rxtr*(T_in-y[2])-DH_fwd/rho/c_p*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        DH_fwd/rho/c_p*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2)+(heat_nom+M_step)/rho/V_rxtr/c_p]
    prob=ODEProblem(f,[CA_nom,T_nom],(0.0,tmax))
    soln=solve(prob,Rosenbrock23())
    t=soln.t
    CA_soln=zeros(size(soln.u)[1])
    T_soln=zeros(size(soln.u)[1])
    for i=1:size(soln.u)[1]
        CA_soln[i]=soln.u[i][1]
        T_soln[i]=soln.u[i][2]
    end
    plot(t,[CA_soln,T_soln],layout=(2,1),legend=false,xlabel="Time (s)", ylabel=["Conc. of A (M)" "Temp. (K)"],gridalpha=0.75,minorgrid=true,minorgridalpha=0.3)
    #    return t,CA_soln,T_soln
end
