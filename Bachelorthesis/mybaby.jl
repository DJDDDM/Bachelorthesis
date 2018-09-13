#!/opt/julia/0.6.4/bin/julia 
using Compat
using TensorToolbox
using TensorOperations
using PyCall
@pyimport f90nml
using TensorOperations

function main()
    two, H, E_nucc, NORB, NELEC = FCIread()
    println("E_nucc= ",E_nucc," NORB= ",NORB," NELEC= ",NELEC)
    NOCC = convert(Int,NELEC/2)
    E_RHF,vals, C = RHF(two,H,E_nucc,NORB,NOCC)
    @time RHF(two,H,E_nucc,NORB,NOCC)

    twoMO = transformtwo(C,two)
    println("transform")
    @time transformtwo(C,two)

    E_MP2 = MP2(twoMO,vals,NOCC,NORB)
    println("E_MP2=",E_MP2)
    @time MP2(two,vals,NOCC,NORB)
end

function FCIread()
    file = open(ARGS[1])
        nml = open("tmpnml.tmp","a")
            while true
                line = readline(file)
                write(nml,line)
                write(nml,"\n")
                if contains(line,"&END") || contains(line,"/")
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
    eig = eigfact(F)
    C = eig[:vectors]
    vals = eig[:values]
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
    twoAO = hosvd(twoAO,eps_abs=100)
    println(ttm(twoAO,C))
    @tensoropt twoMO[i,j,k,l] := twoAO[my,ny,la,si]*C'[my,i]*C'[ny,j]*C[la,k]*C[si,l]
    return twoMO
end


function MP2(two,vals,NOCC,NORB)
    E = 0
    for a = 1:NOCC, b = 1:NOCC, r = NOCC+1:NORB, s = NOCC+1:NORB
    E += ((2*(two[a,r,b,s])*two[r,a,s,b])-(two[a,r,b,s]*two[r,b,s,a]))/(vals[a]+vals[b]-vals[r]-vals[s])
end
return E
end

#function DCI()
#    E_DCI = 0
#    while true
#        B =  
#        for r = NOCC+1:NORB, a = 1:NOCC, s = NOCC+1:NORB, b = 1:NOCC
#            B[r,a,s,b] = 2*two[r,a,s,b]-two[r,b,s,a]
#        end
#        D = 
#        M =
#        E_DCI = -B'inv(M)B
#    end
#end
#
#function Det_func()
#   Det = [1:NOCC]
#end
#   
#function compare(Det1,Det2)
#    for i = 1,NOCC
#        if Det1[i]
    







main()
@time main()
