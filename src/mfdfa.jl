"""
# MFDFA

```
mfdfa(timeseries, q, scales, m; plot = true)
```

Julia implementation of the Multifractal Detrended Fluctuation Analysis (MFDFA) algorithm.

## Input arguments:

* `timeseries`: the timeseries to be analyzed
* `scales`: the scales to be used in the analysis
* `qValues`: the q values to be used in the analysis
* `m`: the order of the polynomial to be used in the local trend estimation
* `plot`: whether to plot the q-order scaling function and mfw (default = true)

## Output:

* `Hq`: the q-order Hurst exponent
* `Hq_r2`: the r-squared value of the q-order Hurst exponent
* `tq`: the q-order mass exponent
* `hq`: the q-order singularity exponent
* `Dq`: the q-order dimension
* `Fq`: the q-order scaling function
* `q`: the q values used in the analysis
* `mfw`: the multifractal width

"""

function mfdfa(timeseries, qValues, scales, m; plot = true)

    X=cumsum(timeseries.-mean(timeseries), dims=1);

    RMS_scales = [Float64[] for i=1:length(scales)]
    qRMS = [Float64[] for i=1:length(qValues),j=1:length(scales)]
    Fq = Array{Float64}(undef, length(qValues), length(scales))
    Hq = Float64[]
    Hq_r2 = Float64[]
    qRegLine = []

    for ns=1:length(scales)
        segments = floor(length(X)/scales[ns])
        for v=1:segments
            Index=Int.([(((v-1)*scales[ns])+1):v*scales[ns];])
            C=Polynomials.PolyCompat.polyfit(Index,X[Index],m);
            pfit = Polynomials.PolyCompat.polyval(C,Index)
            push!(RMS_scales[ns],sqrt(mean((X[Index]-pfit).^2)))
        end
        for nq=1:length(qValues)
            qRMS[nq,ns] = RMS_scales[ns].^qValues[nq]
            Fq[nq,ns]=mean(qRMS[nq,ns]) .^(1/qValues[nq])
        end
        Fq[findall(x->x==0,qValues),ns].=exp(0.5*mean(log.(RMS_scales[ns].^2)));
    end
    for nq=1:length(qValues)
        # using GLM: benefit we get R-squared
        C = GLM.lm(@formula(y~x),DataFrame(y=log2.(Fq[nq,:]),x=log2.(scales)))
        push!(Hq,GLM.coef(C)[2])
        push!(qRegLine,GLM.predict(C))
        push!(Hq_r2, GLM.r²(C))
    end

    tq = (Hq.*qValues).-1
    hq = diff(tq)./(qValues[2]-qValues[1])
    Dq = qValues[1:(length(qValues)-1)].*hq-tq[1:(length(tq)-1)]

    mfw = maximum(hq) - minimum(hq);
    
    
    if plot      

        # Create the plot with the data and the curve
        mf_spectrum_plot = Plots.plot(hq, Dq, legend=false, xlabel = "hq", ylabel = "Dq", title="mfw: $(round(mfw, digits=3))", marker = :circle)
        Hq_q_plot = Plots.plot(qValues, Hq, legend=false, xlabel = "q", ylabel = "Hq", marker = :circle)
                
        figure = Plots.plot(Hq_q_plot, mf_spectrum_plot, size = (1000,500), margin = 5Plots.mm)
    else
        figure = 0
    end

    output = Dict("Hq" => Hq,
                  "Hq_r2" => Hq_r2,
                  "tq" => tq,
                  "hq" => hq,
                  "Dq" => Dq,
                  "Fq" => Fq,
                  "q" => qValues,
                  "mfw" => mfw,
                  "figure" => figure)

    return output
end