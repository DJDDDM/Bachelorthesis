#!/apps/JULIA/1.0.1/bin/julia
module HF
using Compat
using TensorOperations
using PyCall
using Combinatorics
using LinearAlgebra
@pyimport f90nml
import Base.print_matrix

mutable struct DIIS_container
diis_error_matrices
diis_fock_matrices
diis_bmat
diis_energy
citers
end

function main()
    global maxdiis

    global doRHF
    global doUHF
    global doCUHF
    global doMP2

    global NORB
    global NELEC
    global NOCC
    global NVIRT
    global E_nucc

    doRHF = 1
    doUHF = 0
    doCUHF = 1
    doMP2 = 0
    maxdiis = 5

    two, H, E_nucc, NORB, NELEC = FCIread()
    NELEC += parse(Int,ARGS[2])
    println("E_nucc= ",E_nucc," NORB= ",NORB," NELEC= ",NELEC)

    if doRHF == 1
        NOCC = convert(Int,NELEC/2)
        NVIRT = NORB-NOCC
        E_RHF,vals, C = RHF(two,H,E_nucc,NORB,NOCC)
        println("finished RHF E = $(E_RHF)")
    end

    if doUHF == 1
        E_UHF, Ca, Cb, valsa, valsb = UHF(H,two)
        println("finished UHF", "E = ", E_UHF)
    end

    if doCUHF == 1
        E_CUHF, Ca, Cb, valsa, valsb = CUHF(H,two)
        println("finished CUHF", "E = ", E_CUHF)
    end
    
    #show(stdout,"text/plain",C)
    if doMP2 == 1
        twoMO = C' * C' * twoAO * C * C
        E_MP2 = MP2(twoMO,vals,NOCC,NORB)
        println("E_MP2=",E_MP2)
        println("E_RHF+E_MP2= ",E_RHF + E_MP2)
    end

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
    C = Matrix{Float64}(I,NORB,NORB)
    iters = 1
    vals = zeros((NORB))

    DIIS_container = DIIS_init()

    while true
    #C = Transformationmatrix
    #P = Densitymatrix
    #H = onelectronmatrix
    #two = twoelectronmatrix
    #Einsteinconvention

    #C -> P
    P = C[:,1:NOCC]*C[:,1:NOCC]'

    #P -> F
    F = H
    F += 2*tensorcontract(two,["my","ny","si","la"],P,["la","si"])
    F -= tensorcontract(two,["my","la","si","ny"],P,["la","si"])

    #E
    E = E_nucc + tensorcontract(P,["my","ny"],(H+F),["ny","my"])[]

    #DIIS
    F,diis_error = DIIS(DIIS_container,F,P,E,iters)
    #println(diis_error)
    println("E=",E,"iters=",iters)

    #F -> C
    eig = LinearAlgebra.eigen(F)
    C = eig.vectors
    vals = eig.values
    idx = sortperm(vals)
    vals = vals[idx]
    C = C[:,idx]


    #loop control
    if  diis_error < 10^-6
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
    Ca = Matrix{Float64}(I,NORB,NORB)
    Cb = Matrix{Float64}(I,NORB,NORB)
    iters = 1
    valsa = zeros((NORB))
    valsb = zeros((NORB))

    DIISa = DIIS_init()
    DIISb = DIIS_init()

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

        #DIIS
        Fa,Fb,diis_error = DIIS(DIISa,(Fa,Fb),(Pa,Pb),E,iters)

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
        if diis_error < 10^-6 
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

function CUHF(H,two)
    DIIS_container = DIIS_init()

    E = 0.0
    E_old = 0.0
    Da = Matrix{Float64}(I,NORB,NORB)
    Db = Matrix{Float64}(I,NORB,NORB)
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
        vals = real.(eig.values)
        idx = sortperm(vals,rev=true)
        U = U[:,idx]

        #F_AO
        Fcs_AO = H + 2*tensorcontract(two,["i","j","k","l"],Pt_AO,["l","k"]) - tensorcontract(two,["i","l","k","j",],Pt_AO,["l","k"])
        #DIIS
        Fcs_AO,diis_error = DIIS(DIIS_container,Fcs_AO,Pt_AO,E,iters)

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
        valsa = real.(eiga.values)
        idx = sortperm(valsa)
        valsa = valsa[idx]
        Ca = Ca[:,idx]
        Da = U * Ca

        eigb = LinearAlgebra.eigen(Fb_cuhf)
        Cb = eigb.vectors
        valsb = real.(eigb.values)
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
        if iters > 3000
            println("not converged after", iters)
            break
        end
        E_old = E
        iters +=1
    end
    return E, Da, Db, valsa, valsb
end

function DIIS_init()
diis_error_matrices = zeros(maxdiis,NORB,NORB)
diis_fock_matrices= zeros(maxdiis,NORB,NORB)
diis_bmat= zeros(maxdiis,maxdiis)
diis_energy= zeros(maxdiis)
citers = 0
container = DIIS_container(diis_error_matrices, diis_fock_matrices, diis_bmat, diis_energy,citers)
return container
end

function EDIIS(container,bsize)
    #rhs
    rhs = zeros(bsize+1)
    rhs[1:bsize] = container.diis_energy[1:bsize]
    rhs[bsize+1] = 1

    #lhs
    lhs = zeros(bsize+1,bsize+1)
    lhs[bsize+1,:] .= 1
    lhs[:,bsize+1] .= 1
    lhs[bsize+1,bsize+1] = 0
    lhs[1:bsize,1:bsize] = container.diis_bmat[1:bsize,1:bsize]

    #solve
    while true
        X = lhs \ rhs
        #remove
        if minimum(X[1:end-1]) < 0
            i = argmin(X[1:end-1])
            lhs = lhs[setdiff(1:end,i),setdiff(1:end,i)]
            rhs = rhs[setdiff(1:end,i)]
            if length(rhs) == 2
                return [1,0]
                break
            end
            continue
        else
            return X
            break
        end
    end
end

function CDIIS(container)
    container.citers += 1
    csize = min(maxdiis,container.citers)

    rhs = zeros(csize+1)
    rhs[csize+1] = -1

    lhs = zeros(csize+1,csize+1)
    lhs[csize+1,:] .= -1
    lhs[:,csize+1] .= -1
    lhs[csize+1,csize+1] = 0
    lhs[1:csize,1:csize] = container.diis_bmat[1:csize,1:csize]

    X = lhs \ rhs
    return X
end

function DIIS(container,F,P,energy,iters)

    bsize = min(iters,maxdiis)

    if doRHF == 1 || doCUHF == 1
        #copy lower
        for k = min(iters,maxdiis):-1:2
            container.diis_error_matrices[k,:,:] = container.diis_error_matrices[k-1,:,:]
            container.diis_fock_matrices[k,:,:] = container.diis_fock_matrices[k-1,:,:]
            container.diis_energy[k] = container.diis_energy[k-1]
        end
        for i = bsize:-1:2, j = bsize:-1:2
            container.diis_bmat[i,j] = container.diis_bmat[i-1,j-1]
        end

        #errormatrix
        error_mat = F*P
        error_mat -= error_mat'

        #store recent
        container.diis_error_matrices[1,:,:] = error_mat
        container.diis_fock_matrices[1,:,:] = F
        container.diis_energy[1] = energy

        #e-vectors
        e_max = maximum(error_mat)
        e_min = minimum(error_mat)
        e_abs = max(abs(e_max),abs(e_min))

    elseif doUHF == 1
        #initialize Fock as tuple of Fa and Fb
        if iters == 1
            container.diis_fock_matrices = (zeros(maxdiis,NORB,NORB),zeros(maxdiis,NORB,NORB))
        end
        #copy lower
        for k = min(iters,maxdiis):-1:2
            container.diis_error_matrices[k,:,:] = container.diis_error_matrices[k-1,:,:]
            container.diis_fock_matrices[1][k,:,:] = container.diis_fock_matrices[1][k-1,:,:]
            container.diis_fock_matrices[2][k,:,:] = container.diis_fock_matrices[2][k-1,:,:]
            container.diis_energy[k] = container.diis_energy[k-1]
        end
        for i = bsize:-1:2, j = bsize:-1:2
            container.diis_bmat[i,j] = container.diis_bmat[i-1,j-1]
        end

        #errormatrix
        error_mat = F[1]*P[1]+F[2]*P[2]
        error_mat -= P[1]*F[1]+P[2]*F[2]

        #store recent
        container.diis_error_matrices[1,:,:] = error_mat
        container.diis_fock_matrices[1][1,:,:] = F[1][:,:]
        container.diis_fock_matrices[2][1,:,:] = F[2][:,:]
        container.diis_energy[1] = energy
    end

    #e-abs
    e_max = maximum(error_mat)
    e_min = minimum(error_mat)
    e_abs = max(abs(e_max),abs(e_min))

    #bmatrix
    for i = 1:bsize
        container.diis_bmat[1,i] = tr(container.diis_error_matrices[1,:,:]*container.diis_error_matrices[i,:,:])
        container.diis_bmat[i,1] = tr(container.diis_error_matrices[i,:,:]*container.diis_error_matrices[1,:,:])
    end


    #which DIIS?
    #NODIIS
    if e_abs < 10^-6
        X = [1,0]
    #EDIIS
    elseif e_abs > 0.01
        X = EDIIS(container,bsize)
    #CDIIS
    elseif  e_abs <= 0.01
        X = CDIIS(container)

    end
    #new Fock
    if doRHF == 1
        diis_fock = zeros(NORB,NORB)
        for i = 1:length(X)-1
            diis_fock += X[i]*container.diis_fock_matrices[i,:,:]
        end
    return diis_fock, e_abs
    elseif doUHF == 1
        diis_Fa= zeros(NORB,NORB)
        diis_Fb= zeros(NORB,NORB)
        for i = 1:length(X)-1
            diis_Fa += X[i]*container.diis_fock_matrices[1][i,:,:]
            diis_Fb += X[i]*container.diis_fock_matrices[2][i,:,:]
        end
    return diis_Fa,diis_Fb, e_abs
    end
end
end

            

        
import HF
HF.main()

