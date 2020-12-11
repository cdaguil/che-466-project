using ControlSystems
using Polynomials
using Plots
using DifferentialEquations

include("ClosedLoopFunctions.jl")
include("OpenLoopFunctions.jl")

loadProcessData()

#ZN_flow_T(-0.79)
#ZN_heat_C(10500000)

#ZN_flow_C(-137)
#ZN_heat_T(315577)

#openLoopStepHeat(1)

PID_flowC_heatT(0.1,25,9.5,-1/0.00146,0,-0.6911,-1/0.9672,0)
 
