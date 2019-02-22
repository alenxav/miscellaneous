import numpy as np

# gen = numeric matrix where rows are individuals, columns are SNPs
# y = numeric vector with phenotypes

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
