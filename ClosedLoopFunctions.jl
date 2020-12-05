#Functions for testing the closed loop set point performance of PID controllers
#Also, closed loop Ziegler-Nichols experiments
#First order filter functionality added for ISE near-optimal controllers
#Last updated by AA, 2 November 2020

using Plots, DifferentialEquations

function ZN_heat_T(Kc)
    f(y,p,t)=[-flow_nom*y[1]/V_rxtr-2*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        2*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2),
        flow_nom/V_rxtr*(T_in-y[2])-DH_fwd/rho/c_p*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        DH_fwd/rho/c_p*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2)+(heat_nom+Kc*(T_nom+5-y[4]))/rho/V_rxtr/c_p,
        y[2]-y[3],y[3]-y[4]]
    prob=ODEProblem(f,[CA_nom,T_nom,T_nom,T_nom],(0.0,100.0))
    soln=solve(prob,Rosenbrock23())
    t=soln.t
    CA_soln=zeros(size(soln.u)[1])
    T_soln=zeros(size(soln.u)[1])
    for i=1:size(soln.u)[1]
        CA_soln[i]=soln.u[i][1]
        T_soln[i]=soln.u[i][2]
    end
    plot(t,[CA_soln,T_soln],layout=(2,1),legend=false,xlabel="Time (s)", ylabel=["Conc. of A (M)" "Temp. (K)"],gridalpha=0.75,minorgrid=true,minorgridalpha=0.3)
end

function ZN_flow_C(Kc)
    f(y,p,t)=[-(flow_nom+Kc*(CA_nom+0.05-y[4]))*y[1]/V_rxtr-2*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        2*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2),
        (flow_nom+Kc*(CA_nom+0.05-y[4]))/V_rxtr*(T_in-y[2])-DH_fwd/rho/c_p*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        DH_fwd/rho/c_p*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2)+heat_nom/rho/V_rxtr/c_p,
        y[1]-y[3],y[3]-y[4]]
    prob=ODEProblem(f,[CA_nom,T_nom,CA_nom,CA_nom],(0.0,100.0))
    soln=solve(prob,Rosenbrock23())
    t=soln.t
    CA_soln=zeros(size(soln.u)[1])
    T_soln=zeros(size(soln.u)[1])
    for i=1:size(soln.u)[1]
        CA_soln[i]=soln.u[i][1]
        T_soln[i]=soln.u[i][2]
    end
    plot(t,[CA_soln,T_soln],layout=(2,1),legend=false,xlabel="Time (s)", ylabel=["Conc. of A (M)" "Temp. (K)"],gridalpha=0.75,minorgrid=true,minorgridalpha=0.3)
end

function ZN_heat_C(Kc)
    f(y,p,t)=[-flow_nom*y[1]/V_rxtr-2*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        2*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2),
        flow_nom/V_rxtr*(T_in-y[2])-DH_fwd/rho/c_p*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        DH_fwd/rho/c_p*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2)+(heat_nom+Kc*(CA_nom+0.05-y[4]))/rho/V_rxtr/c_p,
        y[1]-y[3],y[3]-y[4]]
    prob=ODEProblem(f,[CA_nom,T_nom,CA_nom,CA_nom],(0.0,1000.0))
    soln=solve(prob,Rosenbrock23())
    t=soln.t
    CA_soln=zeros(size(soln.u)[1])
    T_soln=zeros(size(soln.u)[1])
    for i=1:size(soln.u)[1]
        CA_soln[i]=soln.u[i][1]
        T_soln[i]=soln.u[i][2]
    end
    plot(t,[CA_soln,T_soln],layout=(2,1),legend=false,xlabel="Time (s)", ylabel=["Conc. of A (M)" "Temp. (K)"],gridalpha=0.75,minorgrid=true,minorgridalpha=0.3)
end

function ZN_flow_T(Kc)
    f(y,p,t)=[-(flow_nom+Kc*(T_nom+5-y[4]))*y[1]/V_rxtr-2*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        2*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2),
        (flow_nom+Kc*(T_nom+5-y[4]))/V_rxtr*(T_in-y[2])-DH_fwd/rho/c_p*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        DH_fwd/rho/c_p*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2)+heat_nom/rho/V_rxtr/c_p,
        y[2]-y[3],y[3]-y[4]]
    prob=ODEProblem(f,[CA_nom,T_nom,T_nom,T_nom],(0.0,100.0))
    soln=solve(prob,Rosenbrock23())
    t=soln.t
    CA_soln=zeros(size(soln.u)[1])
    T_soln=zeros(size(soln.u)[1])
    for i=1:size(soln.u)[1]
        CA_soln[i]=soln.u[i][1]
        T_soln[i]=soln.u[i][2]
    end
    plot(t,[CA_soln,T_soln],layout=(2,1),legend=false,xlabel="Time (s)", ylabel=["Conc. of A (M)" "Temp. (K)"],gridalpha=0.75,minorgrid=true,minorgridalpha=0.3)
end

function PID_flowC_heatT(step_C,step_T,Kc_flowC,ti_flowC,td_flowC,Kc_heatT,ti_heatT,td_heatT)
    f(y,p,t)=[-min(1,max(0,(flow_nom+Kc_flowC*(CA_nom+step_C-y[4]+y[7]/ti_flowC-td_flowC*(y[3]-y[4])))))*y[1]/V_rxtr-2*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        2*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2),
        min(1,max(0,(flow_nom+Kc_flowC*(CA_nom+step_C-y[4]+y[7]/ti_flowC-td_flowC*(y[3]-y[4])))))/V_rxtr*(T_in-y[2])-DH_fwd/rho/c_p*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        DH_fwd/rho/c_p*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2)+min(5e5,max(0,(heat_nom+Kc_heatT*(T_nom+step_T-y[6]+y[8]/ti_heatT-td_heatT*(y[5]-y[6])))))/rho/V_rxtr/c_p,
        y[1]-y[3],y[3]-y[4],y[2]-y[5],y[5]-y[6],CA_nom+step_C-y[4],T_nom+step_T-y[6]]
    prob=ODEProblem(f,[CA_nom,T_nom,CA_nom,CA_nom,T_nom,T_nom,0,0],(0.0,10000.0))
    soln=solve(prob,Rosenbrock23())
    t=soln.t
    CA_soln=zeros(size(soln.u)[1])
    T_soln=zeros(size(soln.u)[1])
    auxsoln=zeros(size(soln.u)[1],6)
    for i=1:size(soln.u)[1]
        CA_soln[i]=soln.u[i][1]
        T_soln[i]=soln.u[i][2]
        auxsoln[i,:]=soln.u[i][3:8]
    end
    flow_soln=min.(1,max.(0,flow_nom .+Kc_flowC*(CA_nom .+step_C .-auxsoln[:,2] .+auxsoln[:,5]/ti_flowC .-td_flowC*(auxsoln[:,1] .-auxsoln[:,2]))))
    flow_soln[1]=flow_nom
    heat_soln=min.(5e5,max.(0,heat_nom .+Kc_heatT*(T_nom .+step_T .-auxsoln[:,4] .+auxsoln[:,6]/ti_heatT .-td_heatT*(auxsoln[:,3] .-auxsoln[:,4]))))
    heat_soln[1]=heat_nom
    plot(t,[CA_soln,T_soln,flow_soln,heat_soln],layout=(2,2),legend=false,xlabel="Time (s)", ylabel=["Conc. of A (M)" "Temp. (K)" "Inlet Flowrate (L/s)" "Heating rate (W)"],gridalpha=0.75,minorgrid=true,minorgridalpha=0.3)
end

function PID_flowC_heatT_filter(step_C,step_T,Kc_flowC,ti_flowC,td_flowC,tf_flowC,Kc_heatT,ti_heatT,td_heatT,tf_heatT)
    f(y,p,t)=[-min(1,max(0,y[9]))*y[1]/V_rxtr-2*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        2*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2),
        min(1,max(0,y[9]))/V_rxtr*(T_in-y[2])-DH_fwd/rho/c_p*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        DH_fwd/rho/c_p*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2)+min(5e5,max(0,y[10]))/rho/V_rxtr/c_p,
        y[1]-y[3],y[3]-y[4],y[2]-y[5],y[5]-y[6],CA_nom+step_C-y[4],T_nom+step_T-y[6],
        ((flow_nom+Kc_flowC*(CA_nom+step_C-y[4]+y[7]/ti_flowC-td_flowC*(y[3]-y[4])))-y[9])/tf_flowC,
        ((heat_nom+Kc_heatT*(T_nom+step_T-y[6]+y[8]/ti_heatT-td_heatT*(y[5]-y[6])))-y[10])/tf_heatT]
    prob=ODEProblem(f,[CA_nom,T_nom,CA_nom,CA_nom,T_nom,T_nom,0,0,flow_nom,heat_nom],(0.0,10000.0))
    soln=solve(prob,Rosenbrock23())
    t=soln.t
    CA_soln=zeros(size(soln.u)[1])
    T_soln=zeros(size(soln.u)[1])
    auxsoln=zeros(size(soln.u)[1],6)
    flow_soln=zeros(size(soln.u)[1])
    heat_soln=zeros(size(soln.u)[1])
    for i=1:size(soln.u)[1]
        CA_soln[i]=soln.u[i][1]
        T_soln[i]=soln.u[i][2]
        auxsoln[i,:]=soln.u[i][3:8]
        flow_soln[i]=min(1,max(0,soln.u[i][9]))
        heat_soln[i]=min(5e5,max(0,soln.u[i][10]))
    end
    flow_soln[1]=flow_nom
    heat_soln[1]=heat_nom
    plot(t,[CA_soln,T_soln,flow_soln,heat_soln],layout=(2,2),legend=false,xlabel="Time (s)", ylabel=["Conc. of A (M)" "Temp. (K)" "Inlet Flowrate (L/s)" "Heating rate (W)"],gridalpha=0.75,minorgrid=true,minorgridalpha=0.3)
end

function PID_flowT_heatC(step_C,step_T,Kc_flowT,ti_flowT,td_flowT,Kc_heatC,ti_heatC,td_heatC)
    f(y,p,t)=[-min(1,max(0,(flow_nom+Kc_flowT*(T_nom+step_T-y[6]+y[8]/ti_flowT-td_flowT*(y[5]-y[6])))))*y[1]/V_rxtr-2*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        2*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2),
        min(1,max(0,(flow_nom+Kc_flowT*(T_nom+step_T-y[6]+y[8]/ti_flowT-td_flowT*(y[5]-y[6])))))/V_rxtr*(T_in-y[2])-DH_fwd/rho/c_p*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        DH_fwd/rho/c_p*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2)+min(5e5,max(0,(heat_nom+Kc_heatC*(CA_nom+step_C-y[4]+y[7]/ti_heatC-td_heatC*(y[3]-y[4])))))/rho/V_rxtr/c_p,
        y[1]-y[3],y[3]-y[4],y[2]-y[5],y[5]-y[6],CA_nom+step_C-y[4],T_nom+step_T-y[6]]
    prob=ODEProblem(f,[CA_nom,T_nom,CA_nom,CA_nom,T_nom,T_nom,0,0],(0.0,10000.0))
    soln=solve(prob,Rosenbrock23())
    t=soln.t
    CA_soln=zeros(size(soln.u)[1])
    T_soln=zeros(size(soln.u)[1])
    auxsoln=zeros(size(soln.u)[1],6)
    for i=1:size(soln.u)[1]
        CA_soln[i]=soln.u[i][1]
        T_soln[i]=soln.u[i][2]
        auxsoln[i,:]=soln.u[i][3:8]
    end
    flow_soln=min.(1,max.(0,flow_nom .+Kc_flowT*(T_nom .+step_T .-auxsoln[:,4] .+auxsoln[:,6]/ti_flowT .-td_flowT*(auxsoln[:,3] .-auxsoln[:,4]))))
    flow_soln[1]=flow_nom
    heat_soln=min.(5e5,max.(0,heat_nom .+Kc_heatC*(CA_nom .+step_C .-auxsoln[:,2] .+auxsoln[:,5]/ti_heatC .-td_heatC*(auxsoln[:,1] .-auxsoln[:,2]))))
    heat_soln[1]=heat_nom
    plot(t,[CA_soln,T_soln,flow_soln,heat_soln],layout=(2,2),legend=false,xlabel="Time (s)", ylabel=["Conc. of A (M)" "Temp. (K)" "Inlet Flowrate (L/s)" "Heating rate (W)"],gridalpha=0.75,minorgrid=true,minorgridalpha=0.3)
end

function PID_flowT_heatC_filter(step_C,step_T,Kc_flowT,ti_flowT,td_flowT,tf_flowT,Kc_heatC,ti_heatC,td_heatC,tf_heatC)
    f(y,p,t)=[-min(1,max(0,y[9]))*y[1]/V_rxtr-2*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        2*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2),
        min(1,max(0,y[9]))/V_rxtr*(T_in-y[2])-DH_fwd/rho/c_p*A_fwd*exp(-EA_fwd/R_gas/y[2])*y[1]^2+
        DH_fwd/rho/c_p*A_rev*exp(-EA_rev/R_gas/y[2])*(CB_in-y[1]/2)+min(5e5,max(0,y[10]))/rho/V_rxtr/c_p,
        y[1]-y[3],y[3]-y[4],y[2]-y[5],y[5]-y[6],CA_nom+step_C-y[4],T_nom+step_T-y[6],
        ((flow_nom+Kc_flowT*(T_nom+step_T-y[6]+y[8]/ti_flowT-td_flowT*(y[5]-y[6])))-y[9])/tf_flowT,
        ((heat_nom+Kc_heatC*(CA_nom+step_C-y[4]+y[7]/ti_heatC-td_heatC*(y[3]-y[4])))-y[10])/tf_heatC]
    prob=ODEProblem(f,[CA_nom,T_nom,CA_nom,CA_nom,T_nom,T_nom,0,0,flow_nom,heat_nom],(0.0,10000.0))
    soln=solve(prob,Rosenbrock23())
    t=soln.t
    CA_soln=zeros(size(soln.u)[1])
    T_soln=zeros(size(soln.u)[1])
    auxsoln=zeros(size(soln.u)[1],6)
    flow_soln=zeros(size(soln.u)[1])
    heat_soln=zeros(size(soln.u)[1])
    for i=1:size(soln.u)[1]
        CA_soln[i]=soln.u[i][1]
        T_soln[i]=soln.u[i][2]
        auxsoln[i,:]=soln.u[i][3:8]
        flow_soln[i]=min(1,max(0,soln.u[i][9]))
        heat_soln[i]=min(5e5,max(0,soln.u[i][10]))
    end
    flow_soln[1]=flow_nom
    heat_soln[1]=heat_nom
    plot(t,[CA_soln,T_soln,flow_soln,heat_soln],layout=(2,2),legend=false,xlabel="Time (s)", ylabel=["Conc. of A (M)" "Temp. (K)" "Inlet Flowrate (L/s)" "Heating rate (W)"],gridalpha=0.75,minorgrid=true,minorgridalpha=0.3)
end
