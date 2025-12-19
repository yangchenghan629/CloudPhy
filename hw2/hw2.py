import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

# parameters
cp=1004 # J K-1 kg-1
Rd=287  # J K-1 kg-1
Rv=461  # J K-1 kg-1
cw=4187 # J K-1 kg-1
Lv=2.5e6 # J kg-1
epsilon=0.622 

# initial states
P0=1015*100 # Pa
T0=32+273.15 # K
Td0=21+273.15 # K

# state1 and state2
P1,P2=910*100,800*100 # Pa


# C-C eaquation
def CCeq(T):
    es0=6.112*100 # Pa
    T0=273.15 # K
    es=es0*np.exp(Lv/Rv*(1/T0-1/T))
    return es #Pa

def formula(T,P,w0):
    es=CCeq(T)
    ws=epsilon*es/(P-es)
    w=min(w0,ws)
    chi=max(0,w0-ws)
    Q=w+chi

    e=w*P/(epsilon+w)
    Pd=P-e

    const=(T/(Pd**(Rd/(cp+Q*cw))))*np.exp((w*Lv/(T*(cp+Q*cw))))
    return const

# for initial state
e0=CCeq(Td0)
w0=epsilon*e0/(P0-e0)
chi0=0
Const=formula(T0,P0,w0)

print(f'State 0\nP={P0/100}  T={T0-273.15:.2f}\nw={w0*1000:.2f}  \u03C7={chi0*1000:.2f}\n')

# for state 1
print('State 1')
T1=scipy.optimize.root_scalar(lambda T:formula(T,P1,w0)-Const,method='secant',x0=T0,x1=T0+10,rtol=1e-6).root

es=CCeq(T1)
ws=epsilon*es/(P1-es)
w=min(w0,ws)
chi=max(0,w0-ws)
print(f'T={T1-273.15:.2f}  w={w*1000:.2f}  \u03C7={chi*1000:.2f}')

# for state 2
print('State 2')
T2=scipy.optimize.root_scalar(lambda T:formula(T,P2,w0)-Const,method='secant',x0=T1,x1=T1+10).root

es=CCeq(T2)
ws=epsilon*es/(P2-es)
w=min(w0,ws)
chi=max(0,w0-ws)
print(f'T={T2-273.15:.2f}  w={w*1000:.2f}  \u03C7={chi*1000:.2f}')




