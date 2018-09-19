#!/opt/julia/1.0.0/bin/julia 
using Compat
using TensorToolbox
using TensorOperations
using PyCall
using Combinatorics
using IterTools
using LinearAlgebra
@pyimport f90nml

function main()
    two, H, E_nucc, NORB, NELEC = FCIread()
    println("E_nucc= ",E_nucc," NORB= ",NORB," NELEC= ",NELEC)
    NOCC = convert(Int,NELEC/2)
    global NORB
    global NOCC
    global NVIRT
    NVIRT = NORB-NOCC
    E_RHF,vals, C = RHF(two,H,E_nucc,NORB,NOCC)
    @time RHF(two,H,E_nucc,NORB,NOCC)

    twoMO = transformtwo(C,two)
    HMO = transformone(C,H)
    println("transform")
    @time transformtwo(C,two)

    E_MP2 = MP2(twoMO,vals,NOCC,NORB)
    println("E_MP2=",E_MP2)
    println("E_RHF+E_MP2= ",E_RHF + E_MP2)
    @time MP2(two,vals,NOCC,NORB)

end

function FCIread()
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

