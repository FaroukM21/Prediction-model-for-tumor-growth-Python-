import numpy as np
import random as r
import matplotlib.pyplot as plt
import numpy.random as npr

def choix_aleatoir(M,TUM,p):#entrée une matrice image (notre grille) et p
##    la probabilité de choix de cellule
    while True:
        for x in range(len(M)):
            for y in range(len(M)):
                if r.random()>=p and  (x,y) in TUM:
                    return(x,y)
def liste_diff(M):#liste des cellules différenciées
    L=[]
    for x in range(len(M)):
        for y in range(len(M)):
            if list(M[x,y])==[0,0,255]:#bleu diff cell
                L.append((x,y))
    return(L)
def liste_dediff(M):#liste des cellules dédifférenciées
    L=[]
    for x in range(len(M)):
        for y in range(len(M)):
            if list(M[x,y])==[255,0,0]:#rouge dediff cell
                L.append((x,y))
    return(L)
def liste_Tcells(M):
    L=[]
    for x in range(len(M)):
        for y in range(len(M)):
            if list(M[x,y])==[0,255,0]: #vert Tcell
                L.append((x,y))
    return(L)
def liste_Dead(M): #liste des cellules mortes
    L=[]
    for x in range(len(M)):
        for y in range(len(M)):
            if list(M[x,y])==[250,250,50]:#jaune dead tumor cell
                L.append((x,y))
    return(L)
def rule_for_invasion(Mat,x,y,TUM):
    V=[]#voisinage de la cellule, voisinage de Von Neumann
    if (x+1,y) not in TUM  :
            V.append((x+1,y))
    if (x-1,y)not in TUM  :
            V.append((x-1,y))
    if (x,y+1) not in TUM  :
            V.append((x,y+1))
    if (x,y-1) not in TUM  :
            V.append((x,y-1))
    if (x+1,y+1) not in TUM  :
            V.append((x+1,y+1))
    if (x-1,y-1) not in TUM  :
            V.append((x-1,y-1))
    if (x+1,y-1) not in TUM  :
            V.append((x+1,y-1))
    if (x-1,y+1) not in TUM  :
            V.append((x-1,y+1))
    return(V)
def rule_for_switch(Mat,x,y,T):
    V=[]#voisinage de la cellule, voisinage de Von Neumann
    if (x+1,y)  in T  :
                        V.append((x+1,y))
    if (x-1,y) in T  :
                        V.append((x-1,y))
    if (x,y+1)  in T  :
                        V.append((x,y+1))
    if (x,y-1)  in T  :
                        V.append((x,y-1))
    if (x+1,y+1)  in T  :
            V.append((x+1,y+1))
    if (x-1,y-1)  in T  :
            V.append((x-1,y-1))
    if (x+1,y-1)  in T  :
            V.append((x+1,y-1))
    if (x-1,y+1)  in T  :
            V.append((x-1,y+1))
    return(V)
    

def cellular_automata(bM,dM,cMM,cMD,bD,dD,cDM,cDD,sMD,sDM,bT,dT,sA,dMT,rdecay,n,p):
    """p: probabilité de choix des cellules proliférantes
       n:nombre des instants avant la simulation"""
    
    Mat=np.zeros((101,101,3),np.uint8)
##MOORE voisinage
    Mat[51,51]=[0,0,255]
    Mat[51,50]=[0,0,255]
    Mat[51,52]=[0,0,255]
    Mat[50,51]=[0,0,255]
    Mat[52,51]=[0,0,255]
    Mat[50,52]=[0,255,0]#une cellule T
    Mat[50,50]=[0,255,0]
    Mat[52,50]=[0,255,0]
##    fig, ax=plt.subplots(4)
    fm=lambda M,D,T : M*(bM+dM+cMM*M+cMD*D+sMD+sA+dMT*T)+sDM*D
    fd=lambda M,D,T : (bD+dD+cDD*D+cDM*M+sDM)*D
    ft=lambda M,D,T : (bT*M+dT)*T
    for i in range(n):
        M=liste_diff(Mat)
        D=liste_dediff(Mat)
        T=liste_Tcells(Mat)
        Dead=liste_Dead(Mat)
        TUM=M+D+T+Dead
        T1=TUM[:]
        m=len(M)
        d=len(D)
        t=len(T)
        S=fm(m,d,t)+fd(m,d,t)+ft(m,d,t)
        while T1!=[] :
            (x,y)=choix_aleatoir(Mat,T1,p)
            T1.remove((x,y))
            r1=npr.rand()
            
            if (x,y)in M:
                if r1<=bM*m/S:
                    V=rule_for_invasion(Mat,x,y,TUM)
                    if V!=[]:
                        (a,b)=r.choice(V)
                        TUM.append((a,b))
                        M.append((a,b))
                        Mat[a,b]=[0,0,255]
                elif r1>=1-sMD*m/S :
                    M.remove((x,y))
                    D.append((x,y))
                    Mat[x,y]=[255,0,0]
                elif r1>=1-(sA*m)/S:

                     V=rule_for_switch(Mat,x,y,T)
                     if V!=[]:#rule for invasion
                        D.append((x,y))
                        M.remove((x,y))
                        Mat[x,y]=[255,0,0]
                elif r1 >=1-(dM+cMM*m+cMD*d)*m/S:
                    M.remove((x,y))
                    Dead.append((x,y))
                    Mat[x,y]=[255,255,50]
                elif r1>=1-dMT*m*t/S:
##              rule for invasion, il faut qu'il y ait un voisin sain pour switch
##il faut qu'il y ait cellule T pour avoir un switch
                    V=rule_for_switch(Mat,x,y,T)
                    if V!=[]:
                        M.remove((x,y))
                        Dead.append((x,y))
                        Mat[x,y]=[250,250,50]

            elif (x,y) in D:
                if r1<=bD*d/S:
##              rule for invasion, il faut qu'il y ait un voisin sain pour se reproduire
                    V=rule_for_invasion(Mat,x,y,TUM)
                    if V!=[]:#rule for invasion
                        (a,b)=r.choice(V)
                        TUM.append((a,b))
                        D.append((a,b))
                        Mat[a,b]=[255,0,0]
                
                elif r1>=1-sDM*d/S:
                        D.remove((x,y))
                        M.append((x,y))
                        Mat[x,y]=[0,0,255]
                elif r1>= 1-(dD+cDD*d+cDM*m)*d/S:
                        D.remove((x,y))
                        Dead.append((x,y))
                        Mat[x,y]=[250,250,50]
            elif (x,y) in T:
                
                if r1<=bT*t*m/S:
                    V=rule_for_invasion(Mat,x,y,TUM)
                    if V!=[]:#rule for invasion
                        (a,b)=r.choice(V)
                        TUM.append((a,b))
                        T.append((a,b))
                        Mat[a,b]=[0,255,0]
                elif r1<= dT*t/S:
                    T.remove((x,y))
                    Dead.append((x,y))
                    Mat[x,y]=[250,250,50]
##            else:
##                
##                if r1<=rdecay:
##                    Dead.remove((x,y))
##                    TUM.remove((x,y))
##                    Mat[x,y]=[0,0,0]
##        if i ==15:
##                ax[0].imshow(Mat)
##        if i ==30:
##                ax[1].imshow(Mat)
##        if i ==45:
##                ax[2].imshow(Mat)
##        if i ==100:
##                ax[3].imshow(Mat)
    plt.imshow(Mat)
    plt.show()
                     




cellular_automata(3,1,0.3,1,3,1,0.1,0.3,0.1,1,8,3,4,8,0.35,301,0.01)








    
