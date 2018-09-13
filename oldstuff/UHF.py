#!/usr/bin/python
#########UHF######
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
  if read_addelec: NELEC+=electroninp() #Number of electrons
  else: NELEC+=add_elec
  if NELEC%2==0: NELECa=NELECb=int(NELEC/2)
  else:
    NELECa=int((NELEC+1)/2)
    NELECb=int((NELEC-1)/2)
  E_0_olda=1
  E_0_oldb=1
  Ga = numpy.zeros((NORB+1, NORB+1))   #Gmatrix
  Gb = numpy.zeros((NORB+1, NORB+1))   #Gmatrix
  Gt = numpy.zeros((NORB+1, NORB+1))   #Gmatrix
  Pa = numpy.zeros((NORB+1, NORB+1))   #Gmatrix
  Pb = numpy.zeros((NORB+1, NORB+1))   #Gmatrix
  Pt = numpy.zeros((NORB+1, NORB+1))   #Gmatrix
  dictionary= dict()

  #Control Loop

          #initial P guess
  for i in range (1, NELECa):
    Pa[i,i]=1
  for i in range (1, NELECb):
    Pb[i,i]=1

  Pt=Pa+Pb    

  while True:
    Ga=gmatrixa(dictionary,NORB,Pa)
    Gb=gmatrixb(dictionary,NORB,Pb)
    Gt=gmatrixt(dictionary,NORB,Pt)
    H=hmatrix(NORB,dictionary)
    Fa=fmatrixa(NORB,Ga,Gt,H)
    Fb=fmatrixb(NORB,Gb,Gt,H)
    valsa,vecsa=eigenvectora(Fa)
    valsb,vecsb=eigenvectorb(Fb)
    Pa=pmatrixa(NORB,NELECa,vecsa,Pa)
    Pb=pmatrixb(NORB,NELECb,vecsb,Pb)
    Pt=pmatrixt(Pa,Pb)
    E_0a=energya(NELECa,valsa)
    #print ('new Energya:',E_0a)
    if abs(E_0a-E_0_olda)<10**-7:
      break
    E_0_olda=E_0a
 
  #Energy calculation

  E=dictionary[0,0,0,0]
  for my in range (1,NORB+1):
    for ny in range (1,my):
      E+=(Pt[ny,my]*H[my,ny]+Pa[ny,my]*Fa[my,ny]+Pb[ny,my]*Fb[my,ny])
    ny=my
    E+=0.5*(Pt[ny,my]*H[my,ny]+Pa[ny,my]*Fa[my,ny]+Pb[ny,my]*Fb[my,ny])

  return (E, 'Energy')

def gmatrixa(dictionary,NORB,Pa):          #alpha-G

  Ga = numpy.zeros((NORB+1, NORB+1))
  
  for my in range (1,NORB+1):
    for ny in range (1, my+1):
      Ga[my,ny]=0
      for la in range (1, NORB+1): 
        for si in range (1,la):
          #print("my,ny,la,si",my,ny,la,si)
          Ga[my,ny]+=(-2)*Pa[la,si]*dictionary[dkey2(my,ny,la,si)] 
        si=la
        #print("my,ny,la,la",my,ny,la,si)
        Ga[my,ny]+=(-1)*Pa[la,si]*dictionary[dkey2(my,ny,la,si)] 
  
  for ny in range (2, NORB+1):    #mirror the offdiagonals
    for my in range (1, ny):
      Ga[my,ny]=Ga[ny,my]
  return Ga;


def gmatrixb(dictionary,NORB,Pb):        #beta-G

  Gb = numpy.zeros((NORB+1, NORB+1))
  
  for my in range (1,NORB+1):
    for ny in range (1, my+1):
      Gb[my,ny]=0
      for la in range (1, NORB+1): 
        for si in range (1,la):
          #print("my,ny,la,si",my,ny,la,si)
          Gb[my,ny]+=(-2)*Pb[la,si]*dictionary[dkey2(my,ny,la,si)] 
        si=la
        #print("my,ny,la,la",my,ny,la,si)
        Gb[my,ny]+=(-1)*Pb[la,si]*dictionary[dkey2(my,ny,la,si)] 
  
  for ny in range (2, NORB+1):    #mirror the offdiagonals
    for my in range (1, ny):
      Gb[my,ny]=Gb[ny,my]
  return Gb;


def gmatrixt(dictionary,NORB,Pt):        #T(alpha+beta)-G

  Gt = numpy.zeros((NORB+1, NORB+1))
  
  for my in range (1,NORB+1):
    for ny in range (1, my+1):
      Gt[my,ny]=0
      for la in range (1, NORB+1): 
        for si in range (1,la):
          #print("my,ny,la,si",my,ny,la,si)
          Gt[my,ny]+=2*Pt[la,si]*dictionary[dkey1(my,ny,la,si)]
        si=la
        #print("my,ny,la,la",my,ny,la,si)
        Gt[my,ny]+=1*Pt[la,si]*dictionary[dkey1(my,ny,la,si)]
  
  for ny in range (2, NORB+1):    #mirror the offdiagonals
    for my in range (1, ny):
      Gt[my,ny]=Gt[ny,my]
  return Gt;

def hmatrix(NORB,dictionary):

  H = numpy.zeros((NORB+1, NORB+1))
  for my in range (1,NORB+1):
    for ny in range (1, my+1):
      H[my,ny]=dictionary[my,ny,0,0]
  
  for ny in range (1, NORB+1):    #mirror the offdiagonals
    for my in range (1, ny):
      H[my,ny]=H[ny,my]
  return H;

def fmatrixa(NORB,Ga,Gt,H):       #alpha-F
  Fa = numpy.zeros((NORB+1, NORB+1))
  for my in range (1,NORB+1):
    for ny in range (1, my+1):
      Fa[my,ny]=H[my,ny]+Ga[my,ny]+Gt[my,ny]
  
  for ny in range (2, NORB+1):    #mirror the offdiagonals
    for my in range (1, ny):
      Fa[my,ny]=Fa[ny,my]
  return Fa;

def fmatrixb(NORB,Gb,Gt,H):      #beta-F
  Fb = numpy.zeros((NORB+1, NORB+1))
  for my in range (1,NORB+1):
    for ny in range (1, my+1):
      Fb[my,ny]=H[my,ny]+Gb[my,ny]+Gt[my,ny]
  
  for ny in range (2, NORB+1):    #mirror the offdiagonals
    for my in range (1, ny):
      Fb[my,ny]=Fb[ny,my]
  return Fb;

def eigenvectora (Fa):           #eigenvectors of F-alpha
  valsa, vecsa = LA.eig(Fa[1:,1:])     #gets rid of 0's and sorts the vecs and vals
  idx=valsa.argsort()[::1]
  valsaa=valsa[idx]
  vecsa=vecsa[:,idx]
  return (valsa,vecsa)

def eigenvectorb (Fb):              #eigenvectors of F-beta
  valsb, vecsb = LA.eig(Fb[1:,1:])     #gets rid of 0's and sorts the vecs and vals
  idx=valsb.argsort()[::1]
  valsb=valsb[idx]
  vecsb=vecsb[:,idx]
  return (valsb,vecsb)

def pmatrixa(NORB,NELECa,vecsa,Pa):     #alpha-P

  for my in range (1, NORB+1):
    for ny in range (1, my+1):
      Pa[my,ny]=0
      for a in range (int(NELECa)):
         Pa[my,ny]+=vecsa[a,my-1]*vecsa[a,ny-1]
  
  for ny in range (2, NORB+1):    #mirror the offdiagonals
    for my in range (1, ny):
      Pa[my,ny]=Pa[ny,my]
  return Pa;

def pmatrixb(NORB,NELECb,vecsb,Pb):       #beta-P

  for my in range (1, NORB+1):
    for ny in range (1, my+1):
      Pb[my,ny]=0
      for a in range (int(NELECb)):
         Pb[my,ny]+=vecsb[a,my-1]*vecsb[a,ny-1]
  
  for ny in range (2, NORB+1):    #mirror the offdiagonals
    for my in range (1, ny):
      Pb[my,ny]=Pb[ny,my]
  return Pb;

def pmatrixt(Pa,Pb):        #T(alpha+beta)-P
  Pt=Pa+Pb
  return Pt

def energya(NELECa,valsa):       #alpha-Energy
  
  E_0a=0
    
  for a in range (int(NELECa)):
    E_0a+=valsa[a]
  return (E_0a)

def energyb(NELECb,valsb):      #beta-Energy
  
  E_0b=0
    
  for a in range (int(NELECb)):
    E_0b+=valsb[a]
  return (E_0b)

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
