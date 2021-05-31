# example script for running ChhabraJensen (direct) method

using DataFrames, CSV, GLM, Plots, Plots.PlotMeasures

include("ChhabraJensen.jl")

dataFile = DataFrame(CSV.File("exampleTS.csv", header = false))

dataSeries = dataFile[:,1] |> vec

# check for negative points
minimum(dataSeries)

#= in this example, we have negative values. ChhabraJensen requires positive values so we need to perform a transform. Per Kelty-Stephen et al.
simply adding by some constant is not sufficient. Per [[citation needed]], we can perform a sigmoid transoform:
=#

dataSeries = 1 ./(1 .+ ℯ .^(1 .* dataSeries))

# parameters
qValues = [-10:10;]
scales = [3:10;]

mf_results = ChhabraJensen(dataSeries, qValues, scales)

# Plots

# mf-spectrum plot
mf_spectrum_plot = plot(mf_results["alpha"], mf_results["falpha"],seriestype=:scatter, legend=false, xlabel = "α(q)", ylabel = "f(q)", title=string("mfw: ", round(mf_results["mfw"];digits=3)))

# q-order plot
Dq_q_plot = plot(mf_results["Dq"], mf_results["q"], seriestype=:scatter, legend=false, xlabel = "q", ylabel = "Dq")

# combined plots
plot(Dq_q_plot, mf_spectrum_plot, size = (1000,500), margin = 25px)
