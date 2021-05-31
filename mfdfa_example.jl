# example script for running ChhabraJensen (direct) method

using DataFrames, CSV, GLM, Plots, Plots.PlotMeasures, Polynomials, Polynomials.PolyCompat

include("mfdfa1.jl")

dataFile = DataFrame(CSV.File("exampleTS.csv", header = false))

dataSeries = dataFile[:,1] |> vec

# check for negative points
minimum(dataSeries)


# parameters
q = LinRange(-5,5,101)
exponents = LinRange(4,10,19)
scale = round.(2 .^exponents)
m = 1

mf_results = MFDFA(dataSeries, scale, q, m)

# Plots

# mf-spectrum plot
mf_spectrum_plot = plot(mf_results["hq"], mf_results["Dq"],seriestype=:scatter, legend=false, xlabel = "hq", ylabel = "Dq", title=string("mfw: ", round(mf_results["mfw"];digits=3)))

# q-order plot
Hq_q_plot = plot(mf_results["q"], mf_results["Hq"], seriestype=:scatter, legend=false, xlabel = "q", ylabel = "Hq")

# combined plots
plot(Hq_q_plot, mf_spectrum_plot, size = (1000,500), margin = 25px)
