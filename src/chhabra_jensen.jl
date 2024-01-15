using DataFrames, CSV, GLM, Unitful

"""
```
ChhabraJensen(timeseries, qValues, scales; plot = true)
````

This function calculates the Chhabra-Jensen multifractal spectrum of a time series.

# Arguments
- `timeseries`: A time series of positive values.
- `qValues`: A range of q-values for the multifractal analysis.
- `scales`: A range of scales for the multifractal analysis.
- `plot`: A boolean value indicating whether or not to plot the results.

# Returns
- A dictionary with the following keys:
    - `alpha`: The Hölder exponents.
    - `Rsqr_alpha`: The R² values for the linear regression of `alpha`.
    - `falpha`: The multifractal spectrum.
    - `Rsqr_falpha`: The R² values for the linear regression of `falpha`.
    - `Dq`: The generalized dimensions.
    - `Rsqr_Dq`: The R² values for the linear regression of `Dq`.
    - `q`: The q-values used for the analysis.
    - `mfw`: The width of the multifractal spectrum.
"""
function ChhabraJensen(timeseries, qValues, scales; plot = true)
    timeseries_trim = length(timeseries) % 2^scales[end]
    timeseries = timeseries[1:(end-timeseries_trim)]

    nq, ns = length(qValues), length(scales)
    Ma, Mf, Md = zeros(nq,ns), zeros(nq,ns), zeros(nq,ns)
    alpha, falpha, Dq = zeros(nq), zeros(nq), zeros(nq)
    Rsqr_alpha, Rsqr_falpha, Rsqr_Dq = zeros(nq), zeros(nq), zeros(nq)

    muScale = -log10.(2 .^ scales)
    timeseriesSummed = sum(timeseries)

    for i in 1:nq
        q = qValues[i]
        for j in 1:ns
            window = 2^scales[j]
            timeseriesReshaped = reshape(timeseries, :, window)
            p = sum(timeseriesReshaped, dims = 1) / timeseriesSummed
            Nor = sum(p.^q)
            Md[i,j] = q <= 1 && q > 0 ? sum(p .* log10.(p)) / Nor : log10(Nor)
            mu = (p.^q) / Nor
            Ma[i,j] = sum(mu .* log10.(p))
            Mf[i,j] = sum(mu .* log10.(mu))
        end

        Ma_mdl = lm(@formula(Ma ~ muScale), DataFrame(Ma = Ma[i,:], muScale = muScale))
        Mf_mdl = lm(@formula(Mf ~ muScale), DataFrame(Mf = Mf[i,:], muScale = muScale))
        Md_mdl = lm(@formula(Md ~ muScale), DataFrame(Md = Md[i,:], muScale = muScale))

        alpha[i], falpha[i], Dq[i] = coef(Ma_mdl)[2], coef(Mf_mdl)[2], coef(Md_mdl)[2] / (q - 1)
        Rsqr_alpha[i], Rsqr_falpha[i], Rsqr_Dq[i] = r2(Ma_mdl), r2(Mf_mdl), r2(Md_mdl)
    end
    
    mfw = maximum(alpha) - minimum(alpha)
    
    if plot
        mf_spectrum_plot = Plots.plot(alpha, falpha, seriestype=:scatter, legend=false, xlabel = "α(q)", ylabel = "f(q)", title="mfw: $(round(mfw, digits=3))")
        Dq_q_plot = Plots.plot(Dq, qValues, seriestype=:scatter, legend=false, xlabel = "q", ylabel = "Dq")
        
        diag_plot = Plots.plot(Dq_q_plot, mf_spectrum_plot, size = (1000,500), margin = 5Plots.mm)
    else
        diag_plot = 0
    end
    

    return Dict("alpha" => alpha,
                "Rsqr_alpha" => Rsqr_alpha,
                "falpha" => falpha,
                "Rsqr_falpha" => Rsqr_falpha,
                "Dq" => Dq,
                "Rsqr_Dq" => Rsqr_Dq,
                "q" => qValues,
                "mfw" => mfw,
                "plot" => diag_plot)
end