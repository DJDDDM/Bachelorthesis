#!/opt/julia/1.0.0/bin/julia
using Combinatorics
using IterTools

NORB = 6
NOCC = 4
NVIRT = NORB - NOCC
global NVIRT
global NORB
global NOCC

A = [1,5,3,4]
A = collect(1:NOCC)
B = [1,2,3,4]
onematrix = rand(1:10,5,5)
Gmatrix = rand(1:10,5,5,5,5)
println(A,B)


function sc_exitationgenerator(Nex)
    occ = collect(combinations(1:NOCC,Nex))
    virt = collect(combinations(NOCC+1 : NORB, Nex))
    excited = collect(IterTools.product(occ,virt))
return excited
end


A = sc_exitationgenerator(2)

function sc_difference(A,B)
    counter = 0
    C = [0 0; 0 0]
    D = [0 0; 0 0]
    prefactor = 0

    i = 1
    j = 1
    for i = 1:length(A[1])
        while true
            if A[1][i] == B[1][j]
                break
            end
            if j == length(B[1])
                counter += 1
                if counter == 3
                    break
                end
                C[1,counter] = A[1][i]
                C[2,counter] = B[1][j]
                D[1,counter] = A[1][i]
                D[2,counter] = B[1][j]
                break
            end
            j += 1
        end
        j = 1
        if counter == 3
            break
        end
    end

    if counter == 3
        return C, counter
    end

    i = 1
    j = 1
    for i = 1:length(A[2])
        while true
            if A[2][i] == B[2][j]
                break
            end
            if j == length(B[2])
                counter += 1
                if counter == 3
                    break
                end
                C[1,counter] = A[2][i]
                C[2,counter] = B[2][j]
                D[1,counter] = A[1][i]
                D[2,counter] = B[1][j]
                break
            end
            j += 1
        end
        j = 1
        if counter == 3
            break
        end
    end

    if counter == 0
        prefactor = 1
    elseif counter == 1
        i = D[1,1] - 1 + D[2,1] - 1
        i = i%2 
        if i == 0
            prefactor = 1
        elseif i == 1
            prefactor = -1
        end
    elseif counter == 2
        i = D[1,1] - 1 + D[2,1] -1 + D[1,2] - 2 + D[2,2] - 2
        i = i%2
        if i == 0
            prefactor = 1
        elseif i == 1
            prefactor = -1
        end
    end
            
    return C, counter, prefactor
end




C, counter, prefactor = sc_difference(A[1],A[2])
println(A[1],A[2])
println(C, counter,"pre=",prefactor)


onematrix = rand(1:10,NORB,NORB)
twomatrix = rand(1:10,NORB,NORB,NORB,NORB)

function sc_oneoperator(counter,C,A,onematrix)

    if counter == 0
        value = 0
        j = 1
        for i = 1:NOCC
            if j <= length(A[1]) && i == A[1][j]
                i = A[2][j]
                j += 1
            end
            value += onematrix[i,i]
        end
    elseif counter == 1
        value = onematrix[C[1,1],C[2,1]]
    else
        value = 0
    end
return value
end

value = sc_oneoperator(counter,C,A[1],onematrix)
println(onematrix)
println(value)


function sc_twooperator(counter,C,A,Gmatrix)

    if counter == 0
        k = 1
        l = 1
        value = 0
        for i = 1:NOCC
            for j = 1:NOCC
                if k <= length(A[1]) && i == A[1][k]
                    i = A[2][k]
                    k += 1
                end
                if l <= length(A[1]) && j == A[1][l]
                    j = A[2][l]
                    l +=1
                end
                value += Gmatrix[i,i,j,j]
            end
            l = 1
        end
        value = 0.5 * value
    elseif counter == 1
        value = 0
        k = 1
        for i = 1:NOCC
            if k <= length(A[1]) && i == A[1][k] 
                i = A[2][k]
                k += 1
            end
            value += Gmatrix[C[1,1],C[2,1],i,i]
        end
    elseif counter == 2
        value = Gmatrix[C[1,1],C[2,1],C[1,2],C[2,2]]
    else
        value = 0
    end
return value
end

function sc_onegen(onematrix,A,B)
    C, Counter, Prefactor = sc_difference(A,B)
    value = Prefactor * sc_oneoperator(Counter,C,A,onematrix)
    return value
end

function sc_twogen(twomatrix,A,B)
    C, Counter, Prefactor = sc_difference(A,B)
    value = Prefactor * sc_twooperator(Counter,C,A,twomatrix)
    return value
end

function sc_Hgen(onematrix,twomatrix,A,B)
    return sc_onegen(onematrix,A,B) + 0.5 * sc_twogen(twomatrix,A,B)
end





value = sc_twooperator(counter,C,A[1],twomatrix)
println(value)
println(prefactor*value)

function CID(onematrix,twomatrix)
E_0 = 1.5
E_old = 0.0
E = 0.0
iters = 1
println("onematrix=",onematrix)

C_old = zeros(NORB,NORB,NORB,NORB)
C_neu = copy(C_old)

while true
    #E_corr
    E_corr = 0
    for d = 1:NOCC, c = 1:d, u = NOCC+1:NORB, t = NOCC+1:u
        E_corr += C_neu[c,d,t,u] * twomatrix[c,d,t,u]
    end
    E_corr = E_corr / 2
    E = E_0 - E_corr
    println("CID Energy = ",E)
    
    #C_neu
    for b = 1:NOCC, a = 1:b, s = NOCC+1:NORB, r = NOCC+1:s
        A = [a b; r s]
        value = 0
        for d = 1:NOCC, c = 1:d, u = NOCC+1:NORB, t = NOCC+1:u
            if d == b && c == a && u == s && t == r
                continue
            end
            B = [c d; t u]
            value += C_old[c,d,t,u] * sc_Hgen(onematrix, twomatrix, A, B)
        end
        C_neu[a,b,r,s] = -1 * (twomatrix[r,s,a,b] + value) / (sc_Hgen(onematrix, twomatrix, A, A) - E)
    end
    
    #Loop-control
    if abs(E-E_old) < 10.0^-7
        println("CID converged after",iters,"iterations with energy", E)
        break
    end
    if iters > 1000
        break
    end
    E_old = E
    iters += 1
    C_old= copy(C_neu)
end
println("E= ",E)
end

CID(onematrix,twomatrix)
