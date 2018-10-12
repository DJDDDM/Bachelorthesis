#!/opt/julia/1.0.0/bin/julia
using Compat
using TensorToolbox
using TensorOperations
using PyCall
using Combinatorics
using IterTools
using LinearAlgebra
@pyimport f90nml
import Base.print_matrix

function main()
    two, H, E_nucc, NORB, NELEC = FCIread()
    NELEC += parse(Int,ARGS[2])
    println("E_nucc= ",E_nucc," NORB= ",NORB," NELEC= ",NELEC)
    global NORB
    global NELEC
    global NOCC
    global NVIRT
    global E_nucc
    #NOCC = convert(Int,NELEC/2)
    #NVIRT = NORB-NOCC
    #E_RHF,vals, C = RHF(two,H,E_nucc,NORB,NOCC)
    #@time RHF(two,H,E_nucc,NORB,NOCC)

    #E_RHF,vals, C = RHF_NO(two,H,E_nucc,NORB,NOCC)


    #E_UHF, Ca, Cb, valsa, valsb = UHF(H,two)
    #println("finished UHF", "E = ", E_UHF)

    #E_CUHF, Ca, Cb, valsa, valsb = CUHF_2(H,two)
    #println("finished CUHF", " E = ", E_CUHF)

    E_UHF, Ca, Cb, valsa, valsb = CUHF_3(H,two)
    println("finished CUHF_3", "E = ", E_UHF)

    #E_UHF, Ca, Cb, valsa, valsb = CUHF(H,two)
    #println("finished CUHF", "E = ", E_UHF)

    #3twoMO = transformtwo(C,two)
    #3HMO = transformone(C,H)
    #3println("transform")
    #3@time transformtwo(C,two)

    #3E_MP2 = MP2(twoMO,vals,NOCC,NORB)
    #3println("E_MP2=",E_MP2)
    #3println("E_RHF+E_MP2= ",E_RHF + E_MP2)
    #3@time MP2(two,vals,NOCC,NORB)

end

function FCIread()
    println(ARGS[1])
    file = open(ARGS[1])
        nml = open("tmpnml.tmp","a")
            while true
                line = readline(file)
                write(nml,line)
                write(nml,"\n")
                if occursin("&END",line) || occursin("/",line)
                    break
                end
            end
        close(nml)

        ##call readnml
        nml = f90nml.read("tmpnml.tmp")
        NORB = nml["fci"]["norb"]
        NELEC = nml["fci"]["nelec"]
        rm("tmpnml.tmp")

        #read the integrals
        two = zeros((NORB,NORB,NORB,NORB))
        H = zeros((NORB,NORB))
        E_nucc = 0.0
        while true
            line = readline(file)
            line = split(line)
            integral = parse(Float64,line[1])
            my = parse(Int,line[2])
            ny = parse(Int,line[3])
            la = parse(Int,line[4])
            si = parse(Int,line[5])
            if my != 0 && ny != 0 && la != 0 && si !=0
                two[my,ny,la,si] = integral
                two[my,ny,si,la] = integral
                two[ny,my,si,la] = integral
                two[ny,my,la,si] = integral
                two[la,si,my,ny] = integral
                two[la,si,ny,my] = integral
                two[si,la,ny,my] = integral
                two[si,la,my,ny] = integral
                continue
            elseif my != 0 && ny !=0 && la == 0 && si ==0
                H[my,ny] = integral
                H[ny,my] = integral
                continue
            elseif my == ny == la == si == 0
                E_nucc = integral
                break
                end
            end 
    close(file)
    return two,H,E_nucc,NORB,NELEC
end

function RHF(two,H,E_nucc,NORB,NOCC)
    E = 0.0
    E_old = 0.0
    C = eye(NORB)
    iters = 1
    vals = zeros((NORB))
    

    while true
    #C -> P
    P = C[:,1:NOCC]*C[:,1:NOCC]'

    #P -> F
    F = H
    F += 2*tensorcontract(two,["my","ny","si","la"],P,["la","si"])
    F -= tensorcontract(two,["my","la","si","ny"],P,["la","si"])

    
    #F -> C
    eig = LinearAlgebra.eigen(F)
    C = eig.vectors
    vals = eig.values
    idx = sortperm(vals)
    vals = vals[idx]
    C = C[:,idx]

    #E
    E = E_nucc + tensorcontract(P,["my","ny"],(H+F),["ny","my"])[]
    println("E=",E,"iters=",iters)

    #loop control
    if abs(E-E_old) < 10.0^-7
        println("RHF converged after",iters,"iterations with energy", E)
        break
        end
    if iters > 1000
        break
    end
    E_old = E
    iters += 1
    end
return E, vals, C
end

function transformtwo(C,twoAO)
    @tensor begin
    twoMO[i,j,k,l] := C[i,my]*C[j,ny]*twoAO[my,ny,la,si]*C[la,k]*C[si,l]
    end
    return twoMO
end

function transformone(C,oneAO)
    @tensor begin
    oneMO[i,j] := C[i,my]*oneAO[my,ny]*C[ny,j]
    end
    return oneMO
end


function MP2(two,vals,NOCC,NORB)
    E = 0
    for a = 1:NOCC, b = 1:NOCC, r = NOCC+1:NORB, s = NOCC+1:NORB
    E += ((2*(two[a,r,b,s])*two[r,a,s,b])-(two[a,r,b,s]*two[r,b,s,a]))/(vals[a]+vals[b]-vals[r]-vals[s])
    end
return E
end

function UHF(H,two)

    E = 0.0
    E_old = 0.0
    Ca = eye(NORB)
    Cb = eye(NORB)
    iters = 1
    valsa = zeros((NORB))
    valsb = zeros((NORB))

    if NELEC%2 == 0
        NELECa = convert(Int,(NELEC/2))
        NELECb = convert(Int,(NELEC/2))
    else
        NELECa = convert(Int,((NELEC+1)/2))
        NELECb = convert(Int,((NELEC-1)/2))
    end

    println("NELECa= ",NELECa," NELCb= ", NELECb)

    while true
        #C->P
        Pa = Ca[:,1:NELECa]*Ca[:,1:NELECa]'
        Pb = Cb[:,1:NELECb]*Cb[:,1:NELECb]'
        Pt = Pa + Pb
        #for i = 1:NORB
        #    println(Pt[i,:])
        #end

        #P -> G
        Ga = tensorcontract(Pa,["la","si"], two, ["my","la","si","ny"])
        Gb = tensorcontract(Pb,["la","si"], two, ["my","la","si","ny"])
        Gt = tensorcontract(Pt,["la","si"], two, ["my","ny","si","la"])
        
        #G -> F
        Fa = H+Gt-Ga
        Fb = H+Gt-Gb
        
        #energy
        E = E_nucc + 0.5 * (tensorcontract(Pt, ["ny","my"], H, ["my", "ny"])
                      + tensorcontract(Pa, ["ny","my"], Fa, ["my", "ny"])
                      + tensorcontract(Pb, ["ny","my"], Fb, ["my", "ny"]))

        #F->C
        eiga = LinearAlgebra.eigen(Fa)
        Ca = eiga.vectors
        valsa = eiga.values
        idx = sortperm(valsa)
        valsa = valsa[idx]
        Ca = Ca[:,idx]

        eigb = LinearAlgebra.eigen(Fb)
        Cb = eigb.vectors
        valsb = eigb.values
        idx = sortperm(valsb)
        valsb = valsb[idx]
        Cb = Cb[:,idx]


        #loopcontrol
        println("iters= ",iters,"Energy= ", E)
        if abs(E-E_old)<10.0^-7
            break
        end
        if iters > 1000
            println("not converged after", iters)
            break
        end
        E_old = E
        iters +=1
    end

    println("finished UHF")
    return E, Ca, Cb, valsa, valsb
end

function RHF_NO(two,H,E_nucc,NORB,NOCC)
    E = 0.0
    E_old = 0.0
    D = eye(NORB)
    iters = 1
    vals = zeros((NORB))
    

    while true
    #D: AO->MO
    #C: NO -> MO
    #U: AO -> MO
    #D -> P
    P_AO = D[:,1:NOCC]*D[:,1:NOCC]'

    #NO
    eig = LinearAlgebra.eigen(P_AO)
    U = eig.vectors
    vals = eig.values
    idx = sortperm(vals)
    U = U[:,idx]

    #P -> F
    F_AO = H
    F_AO += 2*tensorcontract(two,["my","ny","si","la"],P_AO,["la","si"])
    F_AO -= tensorcontract(two,["my","la","si","ny"],P_AO,["la","si"])

    #F_NO
    F_NO = U'*F_AO*U

    #F_NO -> C
    eig = LinearAlgebra.eigen(F_NO)
    C = eig.vectors
    vals = eig.values
    idx = sortperm(vals)
    vals = vals[idx]
    C = C[:,idx]

    #D
    D = U*C

    #E
    E = E_nucc + tensorcontract(P_AO,["my","ny"],(H+F_AO),["ny","my"])[]
    println("E=",E,"iters=",iters)

    #loop control
    if abs(E-E_old) < 10.0^-7
        println("RHF_NO converged after",iters,"iterations with energy", E)
        break
        end
    if iters > 1000
        break
    end
    E_old = E
    iters += 1
    end
return E, vals, D
end

function UHF_NO(H,two)

    E = 0.0
    E_old = 0.0
    Da = eye(NORB)
    Db = eye(NORB)
    iters = 1
    valsa = zeros((NORB))
    valsb = zeros((NORB))

    if NELEC%2 == 0
        NELECa = convert(Int,(NELEC/2))
        NELECb = convert(Int,(NELEC/2))
    else
        NELECa = convert(Int,((NELEC+1)/2))
        NELECb = convert(Int,((NELEC-1)/2))
    end

    println("NELECa= ",NELECa," NELCb= ", NELECb)

    while true
        #D->P
        Pa_AO = Da[:,1:NELECa]*Da[:,1:NELECa]'
        Pb_AO = Db[:,1:NELECb]*Db[:,1:NELECb]'
        Pt_AO = Pa_AO + Pb_AO

        #NO
        eig = LinearAlgebra.eigen(Pt_AO)
        U = eig.vectors
        vals = eig.values
        idx = sortperm(vals)
        U = U[:,idx]

        Pa_NO = U' * Pa_AO * U
        Pb_NO = U' * Pa_AO * U
        Pt_NO = Pa_NO + Pb_NO
        
        #P -> G
        Ga = tensorcontract(Pa_AO,["la","si"], two, ["my","la","si","ny"])
        Gb = tensorcontract(Pb_AO,["la","si"], two, ["my","la","si","ny"])
        Gt = tensorcontract(Pt_AO,["la","si"], two, ["my","ny","si","la"])
        
        #G -> F
        Fa_AO = H+Gt-Ga
        Fb_AO = H+Gt-Gb

        #F_NO
        Fa_NO = U' * Fa_AO * U
        Fb_NO = U' * Fb_AO * U
        

        #F->C
        eiga = LinearAlgebra.eigen(Fa_NO)
        Ca = eiga.vectors
        valsa = eiga.values
        idx = sortperm(valsa)
        valsa = valsa[idx]
        Ca = Ca[:,idx]
        Da = U * Ca

        eigb = LinearAlgebra.eigen(Fb_NO)
        Cb = eigb.vectors
        valsb = eigb.values
        idx = sortperm(valsb)
        valsb = valsb[idx]
        Cb = Cb[:,idx]
        Db = U * Cb

        #energy
        E = E_nucc + 0.5 * (tensorcontract(Pt_AO, ["ny","my"], H, ["my", "ny"])
                      + tensorcontract(Pa_AO, ["ny","my"], Fa_AO, ["my", "ny"])
                      + tensorcontract(Pb_AO, ["ny","my"], Fb_AO, ["my", "ny"]))

        #loopcontrol
        println("iters= ",iters,"Energy= ", E)
        if abs(E-E_old)<10.0^-7
            break
        end
        if iters > 100
            println("not converged after", iters)
            break
        end
        E_old = E
        iters +=1
    end
    return E, Da, Db, valsa, valsb
end

function UHF_PM(H,two)

    E = 0.0
    E_old = 0.0
    Da = eye(NORB)
    Db = eye(NORB)
    iters = 1
    valsa = zeros((NORB))
    valsb = zeros((NORB))

    if NELEC%2 == 0
        NELECa = convert(Int,(NELEC/2))
        NELECb = convert(Int,(NELEC/2))
    else
        NELECa = convert(Int,((NELEC+1)/2))
        NELECb = convert(Int,((NELEC-1)/2))
    end

    println("NELECa= ",NELECa," NELCb= ", NELECb)

    while true
        #D->P
        Pa_AO = Da[:,1:NELECa]*Da[:,1:NELECa]'
        Pb_AO = Db[:,1:NELECb]*Db[:,1:NELECb]'
        Pt_AO = (Pa_AO + Pb_AO)/2
        Pd_AO = (Pa_AO - Pb_AO)/2

        #NO
        eig = LinearAlgebra.eigen(Pt_AO)
        U = eig.vectors
        vals = eig.values
        idx = sortperm(vals)
        U = U[:,idx]

        #F_AO
        Fcs_AO = H + 2*tensorcontract(two,["i","j","k","l"],Pt_AO,["l","k"]) - tensorcontract(two,["i","l","k","j",],Pt_AO,["l","k"])
        Duhf_AO = tensorcontract(two,["i","l","k","j"],Pd_AO,["k","l"])

        Fa_AO = Fcs_AO - Duhf_AO
        Fb_AO = Fcs_AO + Duhf_AO
        
        #F_NO
        Fa_NO = U' * Fa_AO * U
        Fb_NO = U' * Fb_AO * U
        

        #F->C
        eiga = LinearAlgebra.eigen(Fa_NO)
        Ca = eiga.vectors
        valsa = eiga.values
        idx = sortperm(valsa)
        valsa = valsa[idx]
        Ca = Ca[:,idx]
        Da = U * Ca

        eigb = LinearAlgebra.eigen(Fb_NO)
        Cb = eigb.vectors
        valsb = eigb.values
        idx = sortperm(valsb)
        valsb = valsb[idx]
        Cb = Cb[:,idx]
        Db = U * Cb

        #energy
        E = (E_nucc
            + 2*tensorcontract(H,["i","j"],Pt_AO,["i","j"])
            + 2*tensorcontract(tensorcontract(two,["i","k","j","l"],Pt_AO,["i","k"]),["j","l"],Pt_AO,["j","l"])[]
            - tensorcontract(tensorcontract(two,["i","l","j","k"],Pt_AO,["i","k"]),["l","j"],Pt_AO,["j","l"])[]
            - tensorcontract(tensorcontract(two,["i","l","j","k"],Pd_AO,["i","k"]),["l","j"],Pd_AO,["j","l"])[])
        
        #loopcontrol
        println("iters= ",iters,"Energy= ", E)
        if abs(E-E_old)<10.0^-7
            break
        end
        if iters > 100
            println("not converged after", iters)
            break
        end
        E_old = E
        iters +=1
    end
    return E, Da, Db, valsa, valsb
end

function CUHF(H,two)

    E = 0.0
    E_old = 0.0
    Da = eye(NORB)
    Db = eye(NORB)
    iters = 1
    valsa = zeros((NORB))
    valsb = zeros((NORB))

    if NELEC%2 == 0
        NELECa = convert(Int,(NELEC/2))
        NELECb = convert(Int,(NELEC/2))
    else
        NELECa = convert(Int,((NELEC+1)/2))
        NELECb = convert(Int,((NELEC-1)/2))
    end

    println("NELECa= ",NELECa," NELCb= ", NELECb)

    while true
        #D->P
        Pa_AO = Da[:,1:NELECa]*Da[:,1:NELECa]'
        Pb_AO = Db[:,1:NELECb]*Db[:,1:NELECb]'
        Pt_AO = (Pa_AO + Pb_AO)/2
        Pd_AO = (Pa_AO - Pb_AO)/2
        #show(stdout, "text/plain",Pt_AO)

        #NO
        eig = LinearAlgebra.eigen(Pt_AO)
        U = eig.vectors
        vals = eig.values
        idx = sortperm(vals)
        U = U[:,idx]
        #show(stdout, "text/plain",U)

        #F_AO
        Fcs_AO = H + 2*tensorcontract(two,["i","j","k","l"],Pt_AO,["l","k"]) - tensorcontract(two,["i","l","k","j",],Pt_AO,["l","k"])
        Duhf_AO = tensorcontract(two,["i","l","k","j"],Pd_AO,["k","l"])

        Fa_AO = Fcs_AO - Duhf_AO
        Fb_AO = Fcs_AO + Duhf_AO

        #F_NO
        Fa_NO = U' * Fa_AO* U
        Fb_NO = U' * Fb_AO* U
        U_cv = zeros(NORB,NORB-(NELECa-NELECb))
        U_cv[:,1:NELECb] = U[:,1:NELECb]
        U_cv[:,NELECb+1:NORB-(NELECa-NELECb)] = U[:,NELECa+1:NORB]
        Duhf_NO = U_cv' * Duhf_AO * U_cv

        Fa_cuhf = copy(Fa_NO)
        Fb_cuhf = copy(Fb_NO)
       
        for i = NELECa:NORB, j = 1:NELECb
            Fa_cuhf[i,j] = Fa_NO[i,j] + Duhf_NO[i-NELECa+1,j]
            Fa_cuhf[j,i] = Fa_NO[j,i] + Duhf_NO[j,i-NELECa+1]

            Fb_cuhf[i,j] = Fb_NO[i,j] - Duhf_NO[i-NELECa+1,j]
            Fb_cuhf[j,i] = Fb_NO[j,i] - Duhf_NO[j,i-NELECa+1]
        end

        ##show(stdout, "text/plain",Fa_cuhf)
        #println("")
        #show(stdout, "text/plain",Fb_AO)
        #println("")
        #show(stdout, "text/plain",Fb_NO)
        #println("")
        #show(stdout, "text/plain",Fb_cuhf)
        
        #F->C
        eiga = LinearAlgebra.eigen(Fa_cuhf)
        Ca = eiga.vectors
        valsa = eiga.values
        idx = sortperm(valsa)
        valsa = valsa[idx]
        Ca = Ca[:,idx]
        Da = U * Ca

        eigb = LinearAlgebra.eigen(Fb_cuhf)
        Cb = eigb.vectors
        valsb = eigb.values
        idx = sortperm(valsb)
        valsb = valsb[idx]
        Cb = Cb[:,idx]
        Db = U * Cb

        #energy
        Edel = 0
            for i = NELECa:NORB, j = 1:NELECb
                Edel += 2 * Duhf_AO[i,j] * Pd_AO[i,j]
                Edel += 2 * Duhf_AO[j,i] * Pd_AO[j,i]
            end
        E = (E_nucc + Edel
            + 2*tensorcontract(H,["i","j"],Pt_AO,["i","j"])
            + 2*tensorcontract(tensorcontract(two,["i","k","j","l"],Pt_AO,["i","k"]),["j","l"],Pt_AO,["j","l"])
            - 1*tensorcontract(tensorcontract(two,["i","l","j","k"],Pt_AO,["i","k"]),["l","j"],Pt_AO,["j","l"])
            - 1*tensorcontract(tensorcontract(two,["i","l","j","k"],Pd_AO,["i","k"]),["l","j"],Pd_AO,["j","l"]))

        
        #loopcontrol
        println("iters= ",iters,"Energy= ", E)
        if abs(E-E_old)<10.0^-7
            break
        end
        if iters > 200
            println("not converged after", iters)
            break
        end
        E_old = E
        iters +=1
        #println("spincontamtination")
        #println(tr(U'*Pa_AO*U*U'*Pb_AO*U))

    end
    return E, Da, Db, valsa, valsb
end

function CUHF_2(H,two)

    E = 0.0
    E_old = 0.0
    Da = eye(NORB)
    Db = eye(NORB)
    iters = 1
    valsa = zeros((NORB))
    valsb = zeros((NORB))

    if NELEC%2 == 0
        NELECa = convert(Int,(NELEC/2))
        NELECb = convert(Int,(NELEC/2))
    else
        NELECa = convert(Int,((NELEC+1)/2))
        NELECb = convert(Int,((NELEC-1)/2))
    end

    println("NELECa= ",NELECa," NELCb= ", NELECb)

    while true
        #D->P
        Pa_AO = Da[:,1:NELECa]*Da[:,1:NELECa]'
        Pb_AO = Db[:,1:NELECb]*Db[:,1:NELECb]'
        Pt_AO = (Pa_AO + Pb_AO)/2
        Pd_AO = (Pa_AO - Pb_AO)/2

        #NO
        eig = LinearAlgebra.eigen(Pt_AO)
        U = eig.vectors
        vals = eig.values
        idx = sortperm(vals)
        U = U[:,idx]

        #F_AO
        Fcs_AO = H + 2*tensorcontract(two,["i","j","k","l"],Pt_AO,["l","k"]) - tensorcontract(two,["i","l","k","j",],Pt_AO,["l","k"])
        Duhf_AO = tensorcontract(two,["i","l","k","j"],Pd_AO,["k","l"])

        Fa_AO = Fcs_AO - Duhf_AO
        Fb_AO = Fcs_AO + Duhf_AO

        #F_NO
        Fa_NO = U' * Fa_AO* U
        Fb_NO = U' * Fb_AO* U
        Fcs_NO = U' * Fcs_AO * U

        Fa_cuhf = copy(Fa_NO)
        Fb_cuhf = copy(Fb_NO)

        #for i = NELECa:NORB, j = 1:NELECb
        #    Fa_cuhf[i,j] = Fcs_NO[i,j]
        #    Fa_cuhf[j,i] = Fcs_NO[j,i]

        #    Fb_cuhf[i,j] = Fcs_NO[i,j]
        #    Fb_cuhf[j,i] = Fcs_NO[j,i]
        #end

        #F->C
        eiga = LinearAlgebra.eigen(Fa_cuhf)
        Ca = eiga.vectors
        valsa = eiga.values
        idx = sortperm(valsa)
        valsa = valsa[idx]
        Ca = Ca[:,idx]
        Da = U * Ca

        eigb = LinearAlgebra.eigen(Fb_cuhf)
        Cb = eigb.vectors
        valsb = eigb.values
        idx = sortperm(valsb)
        valsb = valsb[idx]
        Cb = Cb[:,idx]
        Db = U * Cb

        #energy
        Edel = 0
            for i = NELECa:NORB, j = 1:NELECb
                Edel += 2 * Duhf_AO[i,j] * Pd_AO[i,j]
                Edel += 2 * Duhf_AO[j,i] * Pd_AO[j,i]
            end
        E = (E_nucc + Edel
            + 2*tensorcontract(H,["i","j"],Pt_AO,["i","j"])
            + 2*tensorcontract(tensorcontract(two,["i","k","j","l"],Pt_AO,["i","k"]),["j","l"],Pt_AO,["j","l"])
            - 1*tensorcontract(tensorcontract(two,["i","l","j","k"],Pt_AO,["i","k"]),["l","j"],Pt_AO,["j","l"])
            - 1*tensorcontract(tensorcontract(two,["i","l","j","k"],Pd_AO,["i","k"]),["l","j"],Pd_AO,["j","l"]))

        
        #loopcontrol
        println("iters= ",iters,"Energy= ", E)
        if abs(E-E_old)<10.0^-7
            break
        end
        if iters > 200
            println("not converged after", iters)
            break
        end
        E_old = E
        iters +=1

    end
    return E, Da, Db, valsa, valsb
end

function CUHF_3(H,two)

    E = 0.0
    E_old = 0.0
    Da = eye(NORB)
    Db = eye(NORB)
    iters = 1
    valsa = zeros((NORB))
    valsb = zeros((NORB))

    if NELEC%2 == 0
        NELECa = convert(Int,(NELEC/2))
        NELECb = convert(Int,(NELEC/2))
    else
        NELECa = convert(Int,((NELEC+1)/2))
        NELECb = convert(Int,((NELEC-1)/2))
    end

    println("NELECa= ",NELECa," NELCb= ", NELECb)

    while true
        #D->P
        Pa_AO = Da[:,1:NELECa]*Da[:,1:NELECa]'
        Pb_AO = Db[:,1:NELECb]*Db[:,1:NELECb]'
        Pt_AO = (Pa_AO + Pb_AO)/2
        Pd_AO = (Pa_AO - Pb_AO)/2

        #NO
        eig = LinearAlgebra.eigen(Pt_AO)
        U = eig.vectors
        vals = eig.values
        idx = sortperm(vals,rev=true)
        U = U[:,idx]

        #println("Pt")
        #show(stdout,"text/plain",U'*Pt_AO*U)
        #println("")
        #println("Pd")
        #show(stdout,"text/plain",U'*Pd_AO*U)
        #println("")

        #F_AO
        Fcs_AO = H + 2*tensorcontract(two,["i","j","k","l"],Pt_AO,["l","k"]) - tensorcontract(two,["i","l","k","j",],Pt_AO,["l","k"])
        Duhf_AO = tensorcontract(two,["i","l","k","j"],Pd_AO,["k","l"])

        Fa_AO = Fcs_AO - Duhf_AO
        Fb_AO = Fcs_AO + Duhf_AO
        
        #F_NO
        Fa_NO = U' * Fa_AO * U
        Fb_NO = U' * Fb_AO * U
        Fcs_NO = U' * Fcs_AO * U
        
        Fa_cuhf = copy(Fa_NO)
        Fb_cuhf = copy(Fb_NO)
        for i = NELECa+1:NORB, j = 1:NELECb
            Fa_cuhf[i,j] = Fcs_NO[i,j]
            Fa_cuhf[j,i] = Fcs_NO[j,i]

            Fb_cuhf[i,j] = Fcs_NO[i,j]
            Fb_cuhf[j,i] = Fcs_NO[j,i]
        end

        #F->C
        eiga = LinearAlgebra.eigen(Fa_cuhf)
        Ca = eiga.vectors
        valsa = eiga.values
        idx = sortperm(valsa)
        valsa = valsa[idx]
        Ca = Ca[:,idx]
        Da = U * Ca

        eigb = LinearAlgebra.eigen(Fb_cuhf)
        Cb = eigb.vectors
        valsb = eigb.values
        idx = sortperm(valsb)
        valsb = valsb[idx]
        Cb = Cb[:,idx]
        Db = U * Cb

        #energy
        E = (E_nucc
            + 2*tensorcontract(H,["i","j"],Pt_AO,["i","j"])
            + 2*tensorcontract(tensorcontract(two,["i","k","j","l"],Pt_AO,["i","k"]),["j","l"],Pt_AO,["j","l"])[]
            - tensorcontract(tensorcontract(two,["i","l","j","k"],Pt_AO,["i","k"]),["l","j"],Pt_AO,["j","l"])[]
            - tensorcontract(tensorcontract(two,["i","l","j","k"],Pd_AO,["i","k"]),["l","j"],Pd_AO,["j","l"])[])
        
        #loopcontrol
        println("iters= ",iters,"Energy= ", E)
        if abs(E-E_old)<10.0^-7
            break
        end
        if iters > 300
            println("not converged after", iters)
            break
        end
        E_old = E
        iters +=1
        println("spincontamtination")
        println(tr(U'*Pa_AO*U*U'*Pb_AO*U))
    end
    return E, Da, Db, valsa, valsb
end
main()
