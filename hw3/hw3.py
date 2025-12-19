import numpy as np
import matplotlib.pyplot as plt

class particale:
    def __init__(self,R0,w0):
        self.E=1      
        self.L=2.5e-3  # kg/m3
        self.rho_w=997 # kg/m3
        self.a=8500 # s-1
        self.b=1

        self.R=R0
        self.z=100 # m
        self.w=w0
        self.t=0
    
    def u(self):
        return self.a*self.R**self.b
    
    def step(self,dt):
        uR=self.u()
        dRdt=self.E*self.L/(4*self.rho_w)*uR
        # dRdz=self.E*self.L/(4*self.rho_w)*uR/(self.w-uR)

        self.R+=dRdt*dt
        # self.R+=dRdz*(self.w-uR)*dt
        self.z+=(self.w-uR)*dt
        self.t+=dt

# simulate particle trajectory
def simulate(R0,w0,dt=0.01):
    p=particale(R0,w0)
    height=[p.z]
    radius=[p.R]
    time=[p.t]

    while p.z>0:
        p.step(dt)
        height.append(p.z)
        radius.append(p.R)
        time.append(p.t)
    return np.array(time),np.array(height),np.array(radius)*1e3 # m -> mm

### A
timeA,ZA,RA=simulate(15e-6,2)
### B
timeB,ZB,RB=simulate(30e-6,2)
### C
timeC,ZC,RC=simulate(15e-6,4)
### D
timeD,ZD,RD=simulate(30e-6,4)

max_height=max(np.max(ZA),np.max(ZB),np.max(ZC),np.max(ZD))
max_duration=max(timeA[-1],timeB[-1],timeC[-1],timeD[-1])
max_radius=max(np.max(RA),np.max(RB),np.max(RC),np.max(RD))

outfile=open('hw3_final_radius.txt','w')
outfile.write(f'Radius at cloud base of A : {RA[-1]:.3f} mm\n')
outfile.write(f'Radius at cloud base of B : {RB[-1]:.3f} mm\n')
outfile.write(f'Radius at cloud base of C : {RC[-1]:.3f} mm\n')
outfile.write(f'Radius at cloud base of D : {RD[-1]:.3f} mm\n')
outfile.close()

plt.figure(figsize=(8,6))
plt.title('Z-T Diagram',fontsize=14)
plt.plot(timeA,ZA,'r-',label='A')
plt.plot(timeB,ZB,'b-',label='B')
plt.plot(timeC,ZC,'g--',label='C')
plt.plot(timeD,ZD,'k--',label='D')
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Height above Cloud Base (m)',fontsize=12)
plt.xlim(0,max_duration+100)
plt.ylim(0,max_height+50)
plt.legend()
plt.grid()
plt.savefig('Z_T.png',dpi=300)

plt.figure(figsize=(8,6))
plt.title('Z-R Diagram',fontsize=14)
plt.plot(RA,ZA,'r-',label='A')
plt.plot(RB,ZB,'b-',label='B')
plt.plot(RC,ZC,'g--',label='C')
plt.plot(RD,ZD,'k--',label='D')
plt.xlabel('Radius (mm)',fontsize=12)
plt.ylabel('Height above Cloud Base (m)',fontsize=12)
plt.ylim(0,max_height+50)
plt.xlim(0,max_radius+0.2)
plt.legend()
plt.grid()
plt.savefig('Z_R.png',dpi=300)
