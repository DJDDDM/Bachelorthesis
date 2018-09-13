#!usr/bin/python3

class diis_class:
    def __init__(self,NORB):
        self.diisiters=0
        self.maxdiis=5
        self.diis_error_matrices= numpy.zeros((self.maxdiis, NORB, NORB))
        self.diis_fock_matrices= numpy.zeros_like(self.diis_error_matrices)
        self.diis_energy=numpy.zeros((self.maxdiis))

    def diis_func(self,fock, dens, energy, iters):
        """ Extrapolate new fock matrix based on input fock matrix
            and previous fock-matrices.
            Arguments:
            fock -- current fock matrix
            Returns:
            (fock, error) -- interpolated fock matrix and diis-error
        """
        diis_fock = numpy.zeros_like(fock)
    
        # copy data down to lower storage
        for k in reversed(range(1, min(iters, self.maxdiis))):
            self.diis_error_matrices[k] = self.diis_error_matrices[k-1][:]
            self.diis_fock_matrices[k] = self.diis_fock_matrices[k-1][:]
            self.diis_energy[k] = self.diis_energy[k-1]
    
        # calculate error matrix
        error_mat  = numpy.dot(fock, dens)
        error_mat -= error_mat.T
    
    
        # put data in storage
        self.diis_error_matrices[0]  = error_mat
        self.diis_fock_matrices[0] = fock[:]
        self.diis_energy[0] = energy
    
        #determine wether you go EDIIS or CDIIS
        e_max=error_mat.max()
        e_min=error_mat.min()
        e_abs=max(abs(e_min),abs(e_max))
        print('e_max',e_max, 'e_min',e_min)
        if e_max < 10**-6 and e_min > -10**-6 : #if really small then dont do DIIS
            print('noDIIS')
            return fock,e_abs
        #EDIIS
        if e_max > 0.01 or e_min < -0.01 : #if big then do EDIIS
            print('EDIIS')
            #bmat
            bsize = min(iters, self.maxdiis)
            bmat = 1.0*numpy.ones((bsize+1,bsize+1))
            bmat[bsize, bsize] = 0
            for b1 in range(bsize):
                for b2 in range(bsize):
                    bmat[b1, b2] = numpy.trace(self.diis_error_matrices[b1].dot(self.diis_error_matrices[b2]))
            #rhs
            rhs = numpy.zeros(bsize+1)
            rhs[bsize] = 1
            for a in range(bsize):
                rhs[a]=self.diis_energy[a]
            while True:
                #solve
                C =  numpy.linalg.solve(bmat, rhs)
                #remove C_i<0
                if C[:-1].min() < 0:
                    i=numpy.argmin(C)
                    bmat=numpy.delete(bmat,i,0)
                    bmat=numpy.delete(bmat,i,1)
                    rhs=numpy.delete(rhs,i,0)
                    if len(rhs)==2: #no C_i<0 found
                        print('found no EDIIS')
                        return fock,e_abs
                    continue 
                break

            # form new interpolated diis fock matrix with non0 C_i
            for i, k in enumerate(C[:-1]):
                diis_fock += k*self.diis_fock_matrices[i]

            return diis_fock,e_abs
        #CDIIS
        else:   #if not big then do CDIIS
            print('CDIIS')
            self.diisiters+=1
            bsize = min(self.diisiters, self.maxdiis)
            bmat = -1.0*numpy.ones((bsize+1,bsize+1))
            rhs = numpy.zeros(bsize+1)
            bmat[bsize, bsize] = 0
            rhs[bsize] = -1
            for b1 in range(bsize):
                for b2 in range(bsize):
                    bmat[b1, b2] = numpy.trace(self.diis_error_matrices[b1].dot(self.diis_error_matrices[b2]))
            C =  numpy.linalg.solve(bmat, rhs)

            # form new interpolated diis fock matrix
            for i, k in enumerate(C[:-1]):
                diis_fock += k*self.diis_fock_matrices[i]

            return diis_fock,e_abs
