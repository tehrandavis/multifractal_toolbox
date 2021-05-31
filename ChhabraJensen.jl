#== ChhabraJensen.jl  by Tehran Davis, 2020

This is an Julia-based translation of the Matlab script written by
Authors: Lucas França, Yujiang Wang, José Miranda
https://lucasfr.github.io/, http://xaphire.de/


# --------------------------------------------------------

#INPUT:

This function assumes that your time series is all positive values (sign transform it if needed)

#OUTPUT:
#alpha and falpha scattered against each other give the multifractal spectrum
#qValues and Dq scattered against each other give the generalised dimensions spectrum
#Rsqr_alpha, Rsqr_falpha, and Rsqr_Dq are the R^2 values for each of the values in alpha, falpha, and Dq respectively. It is recommended to run the following code outside of the function:


==#


using DataFrames, CSV, GLM

#timeseries = CSV.read("multifractal_toolbox/exampleTS.txt",header=false)[:,1]
#qValues = [-5:1:5;]
#scales = [3:1:11;]

function ChhabraJensen(timeseries, qValues, scales)

    timeseries_trim = length(timeseries) % 2^scales[length(scales)]
    timeseries = timeseries[1:(length(timeseries)-timeseries_trim)]

    # start function here

    nq=length(qValues)
    ns=length(scales)
    Ma=zeros(nq,ns)
    Mf=zeros(nq,ns)
    Md=zeros(nq,ns)

    muScale= -log10.(2 .^scales)

    alpha=zeros(nq)
    falpha=zeros(nq)
    Dq=zeros(nq)

    Rsqr_alpha=zeros(nq)
    Rsqr_falpha=zeros(nq)
    Rsqr_Dq=zeros(nq)



    for i=1:nq

        q=qValues[i]

        for j=1:ns

            #determine how many windows we will have at this scale
            window=2^scales[j]

            #break the time series into windows & sum
            timeseriesReshaped=reshape(timeseries,:,window)
            timeseriesSummed=sum(timeseries)

            #calculate p
            ps=sum(timeseriesReshaped, dims = 1)
            p=ps./timeseriesSummed


            Nor=sum(p.^q)


            #calculation of Md
            #not accounting for q between 0 and 1
            if q<=1 && q>0
                Md[i,j]=sum(p.*log10.(p))/Nor
            else
                Md[i,j]=log10(Nor)
            end

            #Ma & Mf
            mu=(p.^q)/Nor
            Ma[i,j]=sum(mu.*log10.(p))
            Mf[i,j]=sum(mu.*log10.(mu))
        end

        Ma_mdl = lm(@formula(Ma ~ muScale), DataFrame(Ma = Ma[i,:], muScale = muScale))
        Mf_mdl = lm(@formula(Mf ~ muScale), DataFrame(Mf = Mf[i,:], muScale = muScale))
        Md_mdl = lm(@formula(Md ~ muScale), DataFrame(Md = Md[i,:], muScale = muScale))

        b_Ma = coef(Ma_mdl)[2]
        b_Mf = coef(Mf_mdl)[2]
        b_Md = coef(Md_mdl)[2]

        alpha[i]=b_Ma
        falpha[i]=b_Mf

        if (q<=1 && q>0)
            Dq[i,]=b_Md
        else
            Dq[i]=b_Md/(q-1)
        end

        Rsqr_alpha[i] = r2(Ma_mdl)
        Rsqr_falpha[i] = r2(Mf_mdl)
        Rsqr_Dq[i] = r2(Md_mdl)
    end #looping through scales

    output = Dict("alpha" => alpha,
                    "Rsqr_alpha" => Rsqr_alpha,
                    "falpha" => falpha,
                    "Rsqr_falpha" => Rsqr_falpha,
                    "Dq" => Dq,
                    "Rsqr_Dq" => Rsqr_Dq,
                    "q" => qValues,
                    "mfw" => maximum(alpha)-minimum(alpha)
                    )

    return output

end
