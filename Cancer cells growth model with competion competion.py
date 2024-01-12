import math as m
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as npr
def naissance_mort_competition(b,d,c,T,k):
    Mk=[1]
    Temps=[0]
    delta_t=npr.exponential(1/((b+d+c*k*Mk[-1])*k*Mk[-1]))
    plus_taux=Mk[-1]*b
    minus_taux=Mk[-1]*(d+c*Mk[-1])
    plus=plus_taux/(plus_taux+minus_taux)
    while Temps[-1]+delta_t<=T:
        Temps.append(Temps[-1]+delta_t)
        jump=npr.rand()#parametre pour la methode monte carlo
        if jump < plus:
            Mk.append(Mk[-1]+1/k)
        else:
            Mk.append(Mk[-1]-1/k)
        if Mk[-1]==0:
            break
        detla_t=npr.exponential(1/((b+d+c*k*Mk[-1])*k*Mk[-1]))
        plus_taux=Mk[-1]*b
        minus_taux=Mk[-1]*(d+c*Mk[-1])
        plus=plus_taux/(plus_taux+minus_taux)
    plt.step(Temps,Mk)

naissance_mort_competition(0.12,0.02,0.005,650,1)
##naissance_mort_competition(0.12,0.02,0.005,450,10)
##naissance_mort_competition(0.12,0.02,0.005,450,100)
##naissance_mort_competition(0.12,0.02,0.005,450,1000)
plt.show()

##import math as m
##import matplotlib.pyplot as plt
##import numpy as np
##import numpy.random as npr
##def naissance_mort_competition(b,d,c,T,M0=1,k):
##    M=[M0]
##    Temps=[0]
##    delta_t=npr.exponential(1/((b+d+c*k*M[-1])*k*M[-1]))
##    plus_taux=M[-1]*b
##    minus_taux=M[-1]*(d+c*M[-1])
##    plus=plus_taux/(plus_taux+minus_taux)
##    while Temps[-1]+delta_t<=T:
##        Temps.append(Temps[-1]+delta_t)
##        jump=npr.rand()#parametre pour la methode monte carlo
##        if jump < plus:
##            M.append(M[-1]+1)
##        else:
##            M.append(M[-1]-1)
##        if M[-1]==0:
##            break
##        detla_t=npr.exponential(1/(b+d+c*M[-1]))
##        plus_taux=M[-1]*b
##        minus_taux=M[-1]*(d+c*M[-1])
##        plus=plus_taux/(plus_taux+minus_taux)
##    plt.step(Temps,M)
##
