import scipy.integrate as scint
import numpy as np
import matplotlib.pyplot as plt
m=lambda M,D,T,A : M*(1.9-0.3*M  - 4*A - 28*T )*M + D
d=lambda M,D,T,A : (1.9 - 0.3*D)*D + 0.1*M + 4*A*M
t=lambda M,D,T,A : (8*M - 3 )*T
a=lambda M,D,T,A :-15*A +28*M*T
def ACT(h,tf,M0,D0,T0,A0):
    temps=np.arange(0,tf,h)
    M=[M0]
    D=[D0]
    A=[A0]
    T=[T0]
    n=len(temps)
    for i in range(1,n):
        M.append(M[-1]+h*m(M[-1],D[-1],T[-1],A[-1]))
        D.append(D[-1]+h*d(M[-1],D[-1],T[-1],A[-1]))
        T.append(T[-1]+h*t(M[-1],D[-1],T[-1],A[-1]))
        A.append(A[-1]+h*a(M[-1],D[-1],T[-1],A[-1]))
    X=np.ndarray.tolist(temps)
    plt.plot(X,M)
    plt.plot(X,D)
    plt.plot(X,T)
    plt.plot(X,A)
    plt.show()

