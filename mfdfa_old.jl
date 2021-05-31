#== MFDFA1.jl  by Tehran Davis, 2020

This is an Julia-based translation of the Matlab script based on EAF Ihlen:
Front. Physiol., 04 June 2012 | https://doi.org/10.3389/fphys.2012.00141

==#

using StatsBase, Polynomials, DataFrames, GLM, Polynomials.PolyCompat


function MFDFAold(signal, scale, q, m)
#== Multifractal detrended fluctuation analysis (MFDFA)
[Hq,tq,hq,Dq,Fq]=MFDFA(signal,scale,q,m);
INPUT PARAMETERS---------------------------------------------------------
signal:       input signal
scale:        vector of scales
q:            q-order that weights the local variations
m:            polynomial order for the detrending
Fig:          1/0 flag for output plot of Fq, Hq, tq and multifractal
           spectrum (i.e. Dq versus hq).
OUTPUT VARIABLES---------------------------------------------------------
Hq:           q-order Hurst exponent
tq:           q-order mass exponent
hq:           q-order singularity exponent
Dq:           q-order dimension
Fq:           q-order scaling function
EXAMPLE------------------------------------------------------------------
load fractaldata
scmin=16;
scmax=1024;
scres=19;
exponents=range(log2(scmin),log2(scmax),length = scres);
scale= round.(2 .^exponents)
q=[-5:5;]
m=1
signal1=multifractal;
signal2=monofractal;
signal3=whitenoise;
[Hq1,tq1,hq1,Dq1,Fq1]=MFDFA1(signal1,scale,q,m,1);
[Hq2,tq2,hq2,Dq2,Fq2]=MFDFA1(signal2,scale,q,m,1);
[Hq3,tq3,hq3,Dq3,Fq3]=MFDFA1(signal3,scale,q,m,1);
--------------------------------------------------------------------==#


    X=cumsum(signal.-mean(signal), dims=1);

    RMS_scale = [Float64[] for i=1:length(scale)]
    qRMS = [Float64[] for i=1:length(q),j=1:length(scale)]
    Fq = Array{Float64}(undef, length(q), length(scale))
    Hq = Float64[]
    Hq_r2 = Float64[]
    qRegLine = []


    # if minimum(size(X))!=1||minimum(size(scale))!=1||minimum(size(q))!=1;
    #     error("Input arguments signal, scale and q must be a vector");
    # end
    # if size(X,2)==1;
    #    X=transpose(X);
    # end
    # if min(scale)<m+1
    #    error('The minimum scale must be larger than trend order m+1')
    # end
    for ns=1:length(scale)
        # segments[ns]=floor(length(X)/scale[ns]);
        segments = floor(length(X)/scale[ns])
        for v=1:segments
            Index=Int.([(((v-1)*scale[ns])+1):v*scale[ns];])
            C=Polynomials.fit(Index,X[Index],m);
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

        # using Polynomials - old method
        # C = polyfit(log2.(scale),log2.(Fq[nq,:]),1)
        # push!(Hq,C[1])
        # push!(qRegLine,polyval(C,log2.(scale)))

        # using GLM: benefit we get R-squared
        C = GLM.lm(@formula(y~x),DataFrame(y=log2.(Fq[nq,:]),x=log2.(scale)))
        push!(Hq,GLM.coef(C)[2])
        push!(qRegLine,GLM.predict(C))
        push!(Hq_r2, GLM.r²(C))

    end

    tq = (Hq.*q).-1
    hq = diff(tq)./(q[2]-q[1])
    Dq = q[1:(length(q)-1)].*hq-tq[1:(length(tq)-1)]


    output = Dict("Hq" => Hq,
                    "Hq_r2" => Hq_r2,
                    "tq" => tq,
                    "hq" => hq,
                    "Dq" => Dq,
                    "Fq" => Fq,
                    "q" => q,
                    )

    return output


end
