import numpy as np
import pandas as pd
import csv

Q = np.genfromtxt('Data.csv',delimiter=',')
gen = Q[1:Q.shape[0],1:Q.shape[1]]
y = Q[1:Q.shape[0],0]

#####################
# FITTING THE MODEL #
#####################

iterations = 200
n = gen.shape[0]
p = gen.shape[1]
xx = np.empty(p)
sx = np.empty(p)

for i in range(0,p):
	xx[i] = ((gen[:,i])**2).sum()
	sx[i] = (gen[:,i].sum()**2)/n

MSx = (xx-sx).sum()/(n-1)
mu = y.mean()
m0 = np.empty(1)
e = y-mu
vy = np.dot(y,e)/(n-1)

b = np.repeat(0.0,p)
vb = np.empty(1)
ve = np.dot(y,e)/(n-1)
lmb = MSx

for i in range(0,iterations):	
	for j in range(0,p):
		b0 = b[j]
		b1 = (np.dot(gen[:,j],e)+xx[j]*b0)/(xx[j]+lmb+0.0001)
		e = e-gen[:,j]*(b1-b0)
		b[j] = b1
	m0 = e.mean()
	mu = mu + m0
	e = e - m0	
	ve = np.dot(y,e)/(n-1)
	vb = (vy-ve)/MSx
	lmb = ve/vb

##################
# SUMMARY OUTPUT #
##################

h2 = 1-ve/vy
print('Her H2 = '+str(round(h2,2)))

hat = mu + np.dot(gen,b)
r = np.corrcoef(y,hat)[0,1]
print('GOF (cor) = '+str(round(r2,2)))

##################
# EXPORT RESULTS #
##################

f = open('Data.csv','rb')
reader = csv.reader(f)
headers = reader.next()[1:(p+1)]
f.close()

Final_a = np.r_[mu,b]
Final_b = np.r_[['Mu'],headers]
dataframe = pd.DataFrame(Final_a,index=Final_b,columns=['Effect'])
dataframe.to_csv('Effects.csv')

dataframe = pd.DataFrame(hat,index=np.arange(n)+1,columns=['Hat'])
dataframe.to_csv('GEBV.csv')
