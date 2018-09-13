#!/opt/julia/0.6.4/bin/julia
using Combinatorics
using Iterators

NORB = 6
NOCC = 4
global NORB
global NOCC





A = [1,5,3,4]
A = collect(1:NELEC)
B = [1,2,3,4]
onematrix = rand(1:10,5,5)
Gmatrix = rand(1:10,5,5,5,5)
println(A,B)

sc_exitationgenerator(A,2)

function sc_exitationgenerator(A,Nex)
    occ = collect(combinations(1:NOCC,Nex)
    virt = collect(combinations(NOCC+1:NORB, Nex)
    excited = collect(Iterators.product(occ,virt))

    for i in excited
        println(excited)
    end
end






#function sc_exitationgenerator(Nex)
#
#    A = collect(1:NORB)
#    B = collect(1:NELEC)
#    C = copy(B)
#    D = copy(C)
#    println(C)
#    for i = 1:NELEC, j = NELEC+1:NORB
#            C[i] = A[j]
#            println(C)
#
#            for k = 1:NELEC
#                if k == i
#                    continue
#                end
#                for l = NELEC+1:NORB
#                    if l == j
#                        continue
#                    end
#                    D[k] = A[l]
#                    println(D)
#                    D = copy(C)
#                end
#            end
#            C = copy(B)
#        end
#    end
#
#
#A = collect(1:NORB)
#B = collect(1:NELEC)
#C = copy(B)
#D = copy(A)
#
#j = NELEC+1
#i = 1
#
#while true
#    if D[j] != 0 && C[i] <= NELEC
#        C[i] = D[j]
#        D[j] = 0
#        counter += 1
#        if counter == Nex
#            println(C[i])
#            C = copy(B)
#            D = copy(A)
#            break
#        end
#    end
#end
#
#
#sc_exitationgenerator(1)

function sc_sort(A)
    counter = 0
    i = 1
    while i <= length(A)
        j = i
        while j > 1 && A[j-1] > A[j]
            A[j], A[j-1] = A[j-1], A[j]
            j -= 1
            counter += 1
        end
        i = i+1
    end
    return A, counter
end

A, a = sc_sort(A)
B, b = sc_sort(B)
if (a + b) % 2 == 0
    prefactor = 1
else
    prefactor = -1
end

println(A,B,prefactor)

function sc_compare(A,B)

    counter = 0
    C = [0 0; 0 0]
    i = 1
    j = 1
    while i <= length(A) && j <= length(B)
        if A[i] == B[j] 
            i += 1
            j += 1
        elseif A[i] < B[j]
            counter += 1
            if counter > 2
                break
            end
            C[1,counter] = A[i]
            i +=1
        elseif A[i] > B[j]
            counter += 1
            if counter > 2
                break
            end
            C[2,counter] = B[j]
            j += 1
        end
    end

    if counter < 3
        for i = 1:counter
            if C[1,i] == 0
                C[1,i] = A[end]
            end
            if C[2,i] == 0
                C[2,i] = B[end]
            end
        end
    end

return counter,C
end

counter,C = sc_compare(A,B)
println(counter,C) 

function sc_oneoperator(counter,C,A,onematrix)

    value = 0
    if counter == 0
        for i in A
            value += onematrix[i,i]
        end
    elseif counter == 1
        value = onematrix[C[1,1],C[2,1]]
    else
        value = 0
    end
return value
end
        
function sc_twooperator(counter,C,A,Gmatrix)

    value = 0
    if counter == 0
        for i in A, j in A
            value += Gmatrix[i,i,j,j]
        end
        value = 0.5 * value
    elseif counter == 1
        for i in A
            value += Gmatrix[C[1,1],C[2,1],i,i]
        end
    elseif counter == 2
        value = Gmatrix[C[1,1],C[2,1],C[1,2],C[2,2]]
    else
        value = 0
    end
return value
end








