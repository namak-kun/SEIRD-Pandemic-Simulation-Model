import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def immuneRec(y,t,N,betaS,pid,gamma):
    S,I,R,D = y
    dSdt = float(-S*I*betaS/N)
    dIdt = float(S*I*betaS/N - I*gamma - pid*I)
    dRdt = float(I*gamma)
    dDdt = float(pid*I)
    return dSdt,dIdt,dRdt,dDdt

def immuneRecNot(y,t,N,betaS,betaR,pid,gamma):
    S,I,R,D = y
    dSdt = float(-S*I*betaS/N)
    dIdt = float(S*I*betaS/N + R*I*betaR/N - I*gamma - pid*I)
    dRdt = float(I*gamma - R*I*betaR/N)
    dDdt = float(pid*I)
    return dSdt,dIdt,dRdt,dDdt

N = 1000
psi = 0.2
m = 5
d = 7
betaS = psi*m
pri = 0.1
betaR = pri*m
gamma = 1/d
pid = 0.05

test_type = "RecImmune"
S0, I0 = 999, 1
R0 =  0
D0 = 0

y0 = [S0,I0,R0,D0]

t = np.linspace(0,50,50)
if test_type == "RecImmune":
    sol = sp.integrate.odeint(immuneRec, y0, t, args=(N,betaS,pid,gamma))
elif test_type == "RecNotImmune":
    sol = sp.integrate.odeint(immuneRecNot, y0, t, args=(N,betaS,betaR,pid,gamma))
#print(sol)

plt.plot(t,sol[:,0],'b',label='Susceptible')
plt.plot(t,sol[:,1],'y',label='Infected')
plt.plot(t,sol[:,2],'g',label='recovered')
plt.plot(t,sol[:,3],'',label='dead')
plt.legend()
plt.xlabel('t')
plt.grid()
plt.show()
