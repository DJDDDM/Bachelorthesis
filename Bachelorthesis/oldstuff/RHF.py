#!/usr/bin/python
##############RHF##########
#librarys
import pdb
import sys
from decimal import *
import numpy
from numpy import linalg as LA

def main(add_elec=0,read_addelec=True): 
    #Datadictionary
  NORB = firstline()[0]   #Number of Basisfunctions
  NELEC= firstline()[1]
  diis_N=10          #Number for koefficients in diis
  if read_addelec: NELEC+=electroninp() #Number of electrons
  else: NELEC+=add_elec
  print NELEC
  if (NELEC % 2 != 0) or (NELEC==0):
    print ('odd/zero number of electrons. 1 Electron added')
    NELEC+=1
  E_0_old=1
  G = numpy.zeros((NORB+1, NORB+1))   #Gmatrix
  P = numpy.zeros((NORB+1, NORB+1))   #Pmatrix
  P_cache=numpy.ones((1000,NORB+1,NORB+1))
  dictionary= dict()
  iters=0           #loopindex

  #Control Loop

          #initial P guess
  for x in range (1, int((NELEC/2)+1)):
    P[x,x]=2
  P_cache[0]=P

  while True:
    G=gmatrix(dictionary,NORB,P)
    H=hmatrix(NORB,dictionary)
    F=fmatrix(NORB,G,H)
    vals,vecs=eigenvector(F)
    P=pmatrix(NORB,NELEC,vecs,P)
    #pdb.set_trace()
    #P_cache[iters+1,:,:]=P
    #if iters>15:
    #  P=diis(F,P_cache,iters,diis_N,NORB)
    E_0=energy(vals)
    print ('new Energy:',E_0,'iters',iters)
    if abs(E_0-E_0_old)<10**-5:
      break
    E_0_old=E_0
    iters+=1
 
  #Energy calculation

  E=dictionary[0,0,0,0]
  for my in range (1,NORB+1):
    for ny in range (1,my):
      E+=2*P[my,ny]*(H[my,ny]+0.5*G[my,ny])
    ny=my
    E+=P[my,ny]*(H[my,ny]+0.5*G[my,ny])


  #print ((2*f-g)+dictionary[0,0,0,0],'f-0.5g')
  return (E,'Energy')
  #print (2*h+g+dictionary[0,0,0,0],'h+0.5*g')

def gmatrix(dictionary,NORB,P):

  G = numpy.zeros((NORB+1, NORB+1))
  
  for my in range (1,NORB+1):
    for ny in range (1, my+1):
      G[my,ny]=0
      for la in range (1, NORB+1): 
        for si in range (1,la):
          #print("my,ny,la,si",my,ny,la,si)
          G[my,ny]+=2*P[la,si]*dictionary[dkey1(my,ny,la,si)]
          G[my,ny]+=(-1)*P[la,si]*dictionary[dkey2(my,ny,la,si)] 
        si=la
        #print("my,ny,la,la",my,ny,la,si)
        G[my,ny]+=1*P[la,si]*dictionary[dkey1(my,ny,la,si)]
        G[my,ny]+=(-0.5)*P[la,si]*dictionary[dkey2(my,ny,la,si)] 
  
  for ny in range (2, NORB+1):    #mirror the offdiagonals
    for my in range (1, ny):
      G[my,ny]=G[ny,my]
  return G;

def hmatrix(NORB,dictionary):

  H = numpy.zeros((NORB+1, NORB+1))
  for my in range (1,NORB+1):
    for ny in range (1, my+1):
      H[my,ny]=dictionary[my,ny,0,0]
  
  for ny in range (1, NORB+1):    #mirror the offdiagonals
    for my in range (1, ny):
      H[my,ny]=H[ny,my]
  return H;

def fmatrix(NORB,G,H):
  F = numpy.zeros((NORB+1, NORB+1))
  for my in range (1,NORB+1):
    for ny in range (1, my+1):
      F[my,ny]=H[my,ny]+G[my,ny]
  
  for ny in range (2, NORB+1):    #mirror the offdiagonals
    for my in range (1, ny):
      F[my,ny]=F[ny,my]
  return F;

def eigenvector (F):    

  vals, vecs = LA.eig(F[1:,1:])     #gets rid of 0's and sorts the vecs and vals
  idx=vals.argsort()[::1]
  vals=vals[idx]
  vecs=vecs[:,idx]
  return (vals,vecs)

def pmatrix(NORB,NELEC,vecs,P):

  for my in range (1, NORB+1):
    for ny in range (1, NORB+1):
      P[my,ny]=0
      for a in range (int(NELEC/2)):
         P[my,ny]+=2*vecs[a,my-1]*vecs[a,ny-1]
  
  return P;

def energy(vals):
  
  E_0=0
  for a in range (len(vals)):
    E_0+=vals[a]
  return (E_0)

def dkey1(my,ny,la,si):             #coulomb
  
  if my>la or (my==la and ny>=si):
      return (my,ny,la,si)
  elif my<la or (my==la and ny<si):
      return (la,si,my,ny)
  else: print ('error1')

def dkey2(my,ny,la,si):         #exchange

  if my>=la:
      if ny>=si:
          return(my, la, ny, si)
      elif ny<si:
          return(my,la,si,ny)
      else: print ('error2')
  elif my<la:
      if ny>=si:
          return (la, my,ny,si)
      elif ny<si:
          return(la, my,si,ny)
      else: print ('error3')
  else: print ('error4')

def firstline():          #read NORB and NELEC from FCIDUMP
  
  with open ("FCIDUMP") as f:
      line=f.readline()
      line=line.strip()
      word= line.split('NORB=')[1]
      NORB=(int(word.split(',')[0]))
      word= line.split('NELEC=')[1]
      NELEC=(int(word.split(',')[0]))
      return (NORB,NELEC)

def electroninp():        #read electronaddinginput
  print ('How much electrons you want to add?')
  electronadd=int(input())
  return electronadd

def dict():       #Dictionary

  dictionary={}
  with open ("FCIDUMP") as f:
     for dummy in range(4):
       next(f)
     for line in f:
       integral = float(line.split()[0])   #I_iiii
       first= [line.split()[1], line.split()[2], line.split()[3], line.split()[4]]    #(I_)ijkl
       indizes=[int(first[0]), int(first[1]),int(first[2]), int(first[3])]
       dictionary[(indizes[0]), (indizes[1]), (indizes[2]), (indizes[3])] = integral
  return dictionary



def gtest(NORB,dictionary):   #tests

  G_t = numpy.zeros((NORB+1, NORB+1))

  for my in range (1,NORB+1):
    for ny in range (1,my+1):
      for la in range (1,NORB+1):
        si=la
        G_t[my,ny]+=dictionary[dkey1(my,ny,la,si)]
        G_t[my,ny]+=(-0.5)*dictionary[dkey2(my,ny,la,si)]
        

  for ny in range (2, NORB+1):    #mirror the offdiagonals
    for my in range (1, ny):
      G_t[my,ny]=G_t[ny,my]

  return (G_t)      

def diis(F,P_cache,iters,diis_N,NORB):
  A=numpy.zeros((diis_N+1,diis_N+1))
  K=numpy.zeros((NORB+1,NORB+1))
  Fa=F[1:,1:]
  e=[0]
  B=[1]
  a=-1
  #pdb.set_trace()
  for i in range (diis_N):
    A[i+1,0]=1
    A[0,i+1]=1
    P=P_cache[iters-diis_N+i,:,:]
    P=P[1:,1:]
    e[i]=0
    a=-1
    for my in range (NORB):
      for ny in range (NORB):
        d=0
        for la in range (NORB):
          d+=Fa[my,la]*P[la,ny]-Fa[la,ny]*P[my,la]
        e[i]+=abs(d)
    B.append(0)
    e.append(0)

  for i in range (diis_N):
    for j in range (diis_N):
      A[i+1,j+1]=e[i]*e[j]
  x=[]
  P=numpy.zeros((NORB+1,NORB+1))
  try:
    x=numpy.linalg.solve(A,B)
    for i in range (diis_N):
      P+=x[i+1]*P_cache[iters-diis_N+i]
    return P
  except:
    diis_N-=1
    P=diis(F,P_cache,iters,diis_N,NORB)
    return P


def hf():                   #tests
  E = Decimal(0.)  #Energy
  with open ("FCIDUMP") as f:
    for dummy in range(4):
      next(f)
    for line in f:
      integral = Decimal(line.split()[0])   #I_iiii
      indizes = [int(idx) for idx in line.split()[1:]]    #(I_)ijkl
      if indizes[0] > 5:	 		#outside of Nocc
        continue
      if indizes[0]==indizes[1]==indizes[2]==indizes[3]==0:	#Enuc
        E=E+integral
      elif indizes[0]==indizes[1] and indizes[2]==indizes[3]==0:	#h
        E=E+integral*2
      elif indizes[0]==indizes[1]==indizes[2]==indizes[3]:	#I_iiii
        E=E+integral
      elif indizes[0]==indizes[1] and indizes[2]==indizes[3]:	#J
        E=E+integral*4
      elif indizes[0]==indizes[2] and indizes[1]==indizes[3]:	#K
        E=E-integral*2
  E_hf=E
  return E_hf

if __name__=='__main__': 
  print (main())
