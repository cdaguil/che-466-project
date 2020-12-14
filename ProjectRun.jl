using Polynomials, Plots, DifferentialEquations

include("ClosedLoopFunctions.jl")
include("OpenLoopFunctions.jl")

loadProcessData()
ZN_heat_T(1)
