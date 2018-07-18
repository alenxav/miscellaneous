# loading numpy

import numpy as np
import pandas as pd
import csv

# loading data

Q = np.genfromtxt('Data.csv',delimiter=',')

# splitting file

gen = Q[1:Q.shape[0],1:Q.shape[1]]
y = Q[1:Q.shape[0],0]

# basic settings

df_prior = 4
R2 = 0.5
iterations = 2000
burnin = 500
thinning = 3

# dimensions

n = gen.shape[0]
p = gen.shape[1]

# sum of squares of x

xx = np.empty(p)
for i in range(0,p): xx[i] = ((gen[:,i])**2).sum()

x2 = xx/p
for i in range(0,p): x2[i] = (gen[:,i]).var()

MSx = x2.mean()

# priors

vy = y.var()
Sb_prior = R2 * (df_prior) / MSx
Se_prior = (1-R2) * (df_prior) * vy

# starting parameters

mu = y.mean()
e = y-mu
b = np.empty(p)
ve = 1.1
vb = 1.1
lmb = 1.1

# empty vectors to allocate posterior

Mu = np.empty(1)
B = np.empty(p)
Vb = np.empty(1)
Ve = np.empty(1)

# loop

b0 = np.empty(1)
xy = np.empty(1)
vx = np.empty(1)
m0 = np.empty(1)
MC = 0
print '|',

for i in range(0,iterations):
  
  # update regression coefficients
  
  for j in range(0,p):
    b0 = b[j]
    xy = np.dot(gen[:,j],e)+xx[j]*b0
    vx = xx[j]+lmb    
    b[j] = np.random.normal(xy/vx,ve/vx)
    e = e-gen[:,j]*(b[j]-b0)
    
  # update intercept
  
  m0 = mu
  mu = np.random.normal(m0+e.mean(),ve/n)
  e = e-(mu-m0)
  
  # update variance components
  
  vb = (Sb_prior+np.dot(b,b))/np.random.chisquare(p+df_prior)
  ve = (Se_prior+np.dot(e,e))/np.random.chisquare(n+df_prior)
  lmb = ve/vb
  
  # progress
  
  if(i%(round(iterations/10))==0): print '.',
  
  # save posterior
  
  if (i>=burnin and i%thinning==0):
	  Mu = Mu+mu
	  B = B+b
	  Vb = Vb+vb
	  Ve = Ve+ve
	  MC = MC+1

print '|'

# final estimated

Mu = Mu/MC
B = B/MC
Vb = Vb/MC
Ve = Ve/MC
Va = xx.mean()*Vb

# fitted values

hat = Mu + np.dot(gen,B)
h2 = round(Va/(Va+Ve),4)
r2 = round(np.corrcoef(y,hat)[0,1]**2,4)

# output

f = open('Data.csv','rb')
reader = csv.reader(f)
headers = reader.next()

h2
r2
