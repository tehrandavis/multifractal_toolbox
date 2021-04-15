# mf_plots.jl

#==
use this as a template for plotting mfw plots
==#
using Plots, CSV

qValues = [-5:5;]
scales = [3:11;]

filename = "/Users/tehrandavis/Desktop/multifractal_toolbox/exampleTS.txt"

dataframe = CSV.read(filename, header = false)

ts = dataframe.Column1

df = CAPtoolbox.Multifractal.ChhabraJensen(ts, qValues, scales)

p1 = Plots.plot(df.alpha, df.falpha)
p2 = Plots.plot(df.Dq, qValues)

Plots.plot(p1,p2, legend = false)
