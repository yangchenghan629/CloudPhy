import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

# gamma distribution function
def n(i,D,n0,lam,j):
    return np.array(n0*(D**i)*np.exp(-1*lam*D**j))

# 1) modified Gamma
outfile=open('radar_factor.txt','w+')
N0=1e4
lam=1
j=1
i_array=np.array([-3,0,3])
D=np.linspace(0.1,15,100)
nDs=[]
# plot
fig,ax=plt.subplots(ncols=3,nrows=1,figsize=(8,6))
plt.suptitle('Modified Gamma Distribution',fontsize=14)
for idx,i in enumerate(i_array):
    nDs.append(n(i,D,N0,lam,j))
    ax[idx].set_title(f'i={i}')
    ax[idx].plot(D,nDs[idx])
    ax[idx].set_xlabel('D [mm]',fontsize=12)
    if idx==0:
        ax[0].set_ylabel('n(D) [ $mm^{-i-1}m^{-3}$ ]',fontsize=12)
    ax[idx].grid()
    ax[idx].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    ax[idx].set_yscale('log')
plt.tight_layout()
plt.savefig('gamma_dist.png',dpi=300)
# radar reflectivity factor
for idx in range(3):
    Z=integrate.trapezoid(nDs[idx]*D**6,D)
    outfile.write(f'Z = {Z:.2e} for i = {i_array[idx]}\n')


### 2) bulkwater parameterization
class cloud:
    def __init__(self,N0,i,lam):
        # parameters
        self.fv=1
        self.coef=1e-2
        self.rho=1
        # variables
        self.n0=N0
        self.D=np.linspace(0.1,20,1000)
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

step=20
total_time=60*step
D,nDs,M2s,M3s=simulate(N0=1e6,i=3,lam=1,T=total_time,dt=1)

fig=plt.figure(figsize=(8,6))
plt.title('n(D)',fontsize=14)
for i in range(len(nDs)//600+1):
    plt.plot(D,nDs[i*600],label=f'{10*i} min')
plt.legend()
plt.grid()
plt.xlabel(r'D [$\mu m$]',fontsize=12)
plt.ylabel(r'n(D) [$\mu m^{-4}m^{-3}$]',fontsize=12)
plt.savefig('nD.png',dpi=300)


fig,ax=plt.subplots(figsize=(8,6),ncols=2)
fig.suptitle('Moments',fontsize=14)
for i in range(len(nDs)//600+1):
    ax[0].plot(np.arange(0,total_time+1,1),M2s)
    ax[0].grid()
    ax[0].set_xlabel('Time [sec]',fontsize=12)
    ax[0].set_title('$M_2$',fontsize=12)

    ax[1].plot(np.arange(0,total_time+1,1),M3s)
    ax[1].grid()
    ax[1].set_xlabel('Time [sec]',fontsize=12)
    ax[1].set_title('$M_3$',fontsize=12)
plt.savefig('moments.png',dpi=300)