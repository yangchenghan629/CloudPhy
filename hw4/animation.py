import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from matplotlib.animation import FuncAnimation as F

# gamma distribution function
def n(i,D,n0,lam,j):
    return np.array(n0*(D**i)*np.exp(-1*lam*D**j))

### 2) bulkwater parameterization
class cloud:
    def __init__(self,N0,i,lam):
        # parameters
        self.fv=1
        self.coef=1e-2
        self.rho=1
        # variables
        self.n0=N0
        self.D=np.linspace(0.1,20,100)
        self.i=i
        self.lam=lam
        self.j=1
        self.nD=n(i=self.i,D=self.D,n0=self.n0,lam=self.lam,j=self.j)
        self.M0=self.moments(ith=0)
        self.M2=self.moments(ith=2)
        self.M3=self.moments(ith=3)
    def moments(self,ith):
        func=self.nD*self.D**ith
        return integrate.trapezoid(func,self.D)
    def condense(self,dt):
        dQdt=2*np.pi*self.coef*integrate.trapezoid(self.D*self.nD,self.D)
        dM3dt=6/(np.pi*self.rho)*dQdt
        dM2dt=8/self.rho*self.fv*self.coef*self.M0
        self.M2+=dM2dt*dt
        self.M3+=dM3dt*dt
        self.update()
    def update(self):
        import math
        q=self.M0*self.M3**2/self.M2**3
        self.i=(6-3*q+np.sqrt(q*(q+8)))/(2*(q-1))
        self.lam=self.M2/self.M3*(self.i+3)
        self.n0=self.M0*self.lam**(self.i+1)/math.gamma(self.i+1)
        self.nD=n(self.i,self.D,self.n0,self.lam,self.j)

def simulate(N0,i,lam,T,dt):
    c=cloud(N0,i,lam)
    nDs=[c.nD]
    M2s=[c.M2]
    M3s=[c.M3]
    time=1
    while time<=T:
        c.condense(dt)
        nDs.append(c.nD)
        M2s.append(c.M2)
        M3s.append(c.M3)
        time+=dt
    return c.D,nDs,M2s,M3s

D,nDs,M2s,M3s=simulate(N0=1e6,i=3,lam=1,T=20*60,dt=1)

fig,ax=plt.subplots(figsize=(8,6))
def update_animation(time):
    ax.clear()
    ax.plot(D,nDs[time])
    ax.grid()
    ax.set_title(f'time={time} sec')
    ax.set_xlim(0,20)
    ax.set_ylim(0,3e6)
    ax.set_xlabel(r'D [$\mu m$]',fontsize=12)
    ax.set_ylabel(r'n(D) [$\mu m^{-4}m^{-3}$]',fontsize=12)
    return ax
plt.suptitle('Condensation Growth',fontsize=14)
ani=F(fig,update_animation,frames=len(nDs))
ani.save('animation.mp4',fps=120)