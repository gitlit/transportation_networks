import numpy as np
import pandas as pd
from scipy import optimize 
from scipy.optimize import minimize_scalar
import scipy.linalg as la
import numpy.linalg as la2
import scipy.integrate as integrate
import csv

city = "SiouxFalls"

def convert_trips(input, output):
    origin = -1
    out = []
    with open(input, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row)>0:
                l = row[0].split()
                if l[0] == 'Origin':
                    origin = l[1]
                    if origin != '1':
                        out.append('\n')
                elif origin != -1:
                    for i,e in enumerate(l):
                        if e != '1' and i%3 == 0:
                            out.append(',')
                        if i%3 == 2:
                            out.append('{}'.format(e[:-1]))

    with open(output, "w") as text_file:
        text_file.write(''.join(out))


def convert_net(input, output, output2):
    flag = False
    nodes_number = 0
    links_number = 0
    out = []
    out2 = []
    with open(input, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) > 0 and flag == False and len(row[0].split()) > 3:
                if row[0].split()[0] + row[0].split()[1] + row[0].split()[2] == '<NUMBEROFNODES>':
                    nodes_number = int(row[0].split()[3])
                if row[0].split()[0] + row[0].split()[1] + row[0].split()[2] == '<NUMBEROFLINKS>':
                    links_number = int(row[0].split()[3])

            if len(row) > 0 and flag == True and len(row[0].split()) > 3:
                a = int(row[0].split()[0])
                b = int(row[0].split()[1])
                fft = float(row[0].split()[4])
                cap = float(row[0].split()[2])
                out2.append('{},{}\n'.format(fft, cap))
                c = np.zeros(nodes_number,dtype=np.int)
                c[int(a)-1] = 1
                c[int(b)-1] = -1
                i = 0
                while i < nodes_number - 1:
                    out.append('{},'.format(c[i]))
                    i = i + 1
                out.append('{}\n'.format(c[nodes_number - 1]))

            if len(row) > 0 and len(row[0].split()) > 0 and flag == False:
                if row[0].split()[0] == '~':
                    flag = True

    with open(output, "w") as text_file:
        text_file.write(''.join(out))
    with open(output2, "w") as text_file:
        text_file.write(''.join(out2))


def nk(input):
    flag = False
    with open(input+"_net.tntp", 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) > 0 and flag == False and len(row[0].split()) > 3:
                if row[0].split()[0] + row[0].split()[1] + row[0].split()[2] == '<NUMBEROFNODES>':
                    nodes_number = int(row[0].split()[3])
                if row[0].split()[0] + row[0].split()[1] + row[0].split()[2] == '<NUMBEROFLINKS>':
                    links_number = int(row[0].split()[3])
            if len(row) > 0 and flag == False and row[0].split()[0] == '~':
                flag = True
    nk = [nodes_number, links_number]
    return nk


def convert(input):
    convert_trips(str(input)+"_trips.tntp", "Q.csv")
    convert_net(str(input)+"_net.tntp", "linknode.csv", "coeff.csv")


def compute_time(t0,fe,f0):
	te = t0*(1+0.15*(fe/f0)**4)
	return te


def compute_psi(gamma,fe,f0,t0,ye):
	psi = 0
	for i in range(len(fe)):
		psi += integrate.quad(lambda z: compute_time(t0[i],z,f0[i]),0,fe[i]+gamma*(ye[i]-fe[i]))[0]
	return psi


def linear_search(fe,f0,t0,ye):
	gamma = minimize_scalar(compute_psi, args=(fe, f0, t0,ye), bounds = (0,1), method = 'Bounded')
	return gamma.x


#1. prepare data
convert(city)

LinkNode = pd.read_csv("linknode.csv", header = None)
LinkNode = LinkNode.as_matrix()

Q = pd.read_csv("Q.csv", header = None)
Q = Q.as_matrix()

k = nk(city)[0] # number of total nodes
n = nk(city)[1] # number of total links
s = (n,k)


coeff = pd.read_csv("coeff.csv", header = None)
coeff = coeff.as_matrix()


#2.looking for the shortest path
t0 = coeff [:,0] # free flow travel time
f0 = coeff[:,1] # capacity for each link

origq = np.sum(Q, axis = 1) # row sums of Q: total flow from origin i
destq = -Q 
s2 = (k,k)
RHT = np.zeros(s2)# each row represents an origin, each column represents the flow on node k (with origin i )
RHT = -Q
np.fill_diagonal(RHT, origq)
#print RHT

c0_0 = np.transpose(t0)
c_0 = np.tile(c0_0,k)

A0 = np.transpose(LinkNode) # Construct block matrix for A_eq
A1 = [A0]*k
A = la.block_diag(*A1)

b0 = np.transpose(RHT)
b = np.ravel(b0, order = 'F')[:,np.newaxis] # construct long b_eq

ybounds = (0, None)
result = optimize.linprog(
	c_0, A_eq = A, b_eq = b, bounds = (ybounds), options = {"disp":True, "maxiter":2000,"bland":True} 
	)

result = np.reshape(result['x'],(k,n))
fe = np.sum(result, axis = 0) # intialization fe
te = compute_time(t0,fe,f0)
#end

step = 0
tenorm = 1000000

iteration = []
psi = []

while (tenorm>7.6): # allow each link has 0.1 diff. in te on average
	print "step ", step
	iteration.append(step)
	te_old = te

	c0 = np.transpose(te)
	c = np.tile(c0,k)
	result = optimize.linprog(
	c, A_eq = A, b_eq = b, bounds = (ybounds), options = {"disp":True, "maxiter":2000,"bland":True} 
	)

	resultx = np.reshape(result['x'],(k,n))
	ye = np.sum(resultx, axis = 0) # yn

	gamma = linear_search(fe,f0,t0,ye)
	print "gamma is", gamma
	fe = (1-gamma)*fe + gamma * ye 
	print "fe is ",fe

	te = compute_time(t0,fe,f0)
	tenorm = la2.norm(te-te_old)
	z = np.dot(np.transpose(fe),te)
	psi.append(z)
	print "te is ", te
	print "norm of te is ", tenorm
	step +=1 

import matplotlib.pyplot as plt
plt.plot(iteration,psi,'ro')
plt.xlabel('iteration number')
plt.ylabel('psi(x)')
plt.show()


np.savetxt("te", te, delimiter = ",")
np.savetxt("fe", fe, delimiter = ",")