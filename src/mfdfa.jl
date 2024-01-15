"""
# MFDFA

```
mfdfa(signal, scale, q, m; plot = true)
```

Julia implementation of the Multifractal Detrended Fluctuation Analysis (MFDFA) algorithm.

## Input arguments:

* `signal`: the signal to be analyzed
* `scale`: the scales to be used in the analysis
* `q`: the q values to be used in the analysis
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

function mfdfa(signal, scale, q, m; plot = true)

    X=cumsum(signal.-mean(signal), dims=1);

    RMS_scale = [Float64[] for i=1:length(scale)]
    qRMS = [Float64[] for i=1:length(q),j=1:length(scale)]
    Fq = Array{Float64}(undef, length(q), length(scale))
    Hq = Float64[]
    Hq_r2 = Float64[]
    qRegLine = []

    for ns=1:length(scale)
        segments = floor(length(X)/scale[ns])
        for v=1:segments
            Index=Int.([(((v-1)*scale[ns])+1):v*scale[ns];])
            C=Polynomials.PolyCompat.polyfit(Index,X[Index],m);
            pfit = Polynomials.PolyCompat.polyval(C,Index)
            push!(RMS_scale[ns],sqrt(mean((X[Index]-pfit).^2)))
        end
        for nq=1:length(q)
            qRMS[nq,ns] = RMS_scale[ns].^q[nq]
            Fq[nq,ns]=mean(qRMS[nq,ns]) .^(1/q[nq])
        end
        Fq[findall(x->x==0,q),ns].=exp(0.5*mean(log.(RMS_scale[ns].^2)));
    end
    for nq=1:length(q)
        # using GLM: benefit we get R-squared
        C = GLM.lm(@formula(y~x),DataFrame(y=log2.(Fq[nq,:]),x=log2.(scale)))
        push!(Hq,GLM.coef(C)[2])
        push!(qRegLine,GLM.predict(C))
        push!(Hq_r2, GLM.r²(C))
    end

    tq = (Hq.*q).-1
    hq = diff(tq)./(q[2]-q[1])
    Dq = q[1:(length(q)-1)].*hq-tq[1:(length(tq)-1)]

    mfw = maximum(hq) - minimum(hq);
    
    if plot
        mf_spectrum_plot = Plots.plot(hq, Dq,seriestype=:scatter, legend=false, xlabel = "hq", ylabel = "Dq", title="mfw: $(round(mfw, digits=3))")
        Hq_q_plot = Plots.plot(q, Hq, seriestype=:scatter, legend=false, xlabel = "q", ylabel = "Hq")
        diag_plot = Plots.plot(Hq_q_plot, mf_spectrum_plot, size = (1000,500), margin = 5Plots.mm)
    else
        diag_plot = 0
    end

    output = Dict("Hq" => Hq,
                  "Hq_r2" => Hq_r2,
                  "tq" => tq,
                  "hq" => hq,
                  "Dq" => Dq,
                  "Fq" => Fq,
                  "q" => q,
                  "mfw" => mfw,
                  "diag_plot" => diag_plot)

    return output
end