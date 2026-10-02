#########PID controller + Uncontrolled system

from numpy import zeros,dot,atleast_1d
import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt
from math import sqrt,floor
import pandas as pd
import matplotlib as mpl
################################################################################


def dampparams(alpha,M,K): #return damping matrix 
    
    beta=0
    w1=sqrt(1-(sqrt(2)/2))*sqrt(K[1,1]/M[1,1]) #first natural freq
    w2=sqrt(1+(sqrt(2)/2))*sqrt(K[1,1]/M[1,1]) #second natural freq
    zeta1 = sqrt((2*sqrt(2)-2)/4)*alpha*sqrt(M[1,1]/K[1,1])+sqrt(2/(2-sqrt(2)))*beta*sqrt(M[1,1]/K[1,1])
    zeta2 = sqrt((2*sqrt(2)+2)/4)*alpha*sqrt(M[1,1]/K[1,1])+sqrt(2/(2+sqrt(2)))*beta*sqrt(M[1,1]/K[1,1])
    phiinv11 = (K[0,0]-w2**2*M[0,0])/(M[0,0]*(w1**2-w2**2))
    phiinv12 = K[0,1]/(M[0,0]*(w1**2-w2**2))
    phiinv21 = -(K[0,0]-w1**2*M[0,0])/(M[0,0]*(w1**2-w2**2))
    phiinv22 = -K[0,1]/(M[0,0]*(w1**2-w2**2))
    phiinv=np.array([[phiinv11,phiinv12],[phiinv21,phiinv22]])
    C1=2*zeta1*w1*(M[0,0]+(M[1,1]*(K[0,0]-w1**2*M[0,0])**2)/K[0,1]**2)
    C2=2*zeta2*w2*(M[0,0]+(M[1,1]*(K[0,0]-w2**2*M[0,0])**2)/K[0,1]**2)
    c = np.matmul(np.matmul(np.transpose(phiinv),np.array([[C1,0],[0,C2]])),phiinv)

    return c

##base acceleration 

def ddxg_el(Tf,dt):
    df = pd.read_excel (r'ELR.xlsx')
    baseacc=np.array(df)
    baseacc = np.insert(baseacc, 0, np.array([0,0.0063]), axis=0)
    NT=floor(Tf/dt)
    t=np.linspace(0,Tf,NT)
    new_acc = np.interp(t, baseacc[:,0], baseacc[:,1])*9.80665
    return new_acc

def ddxg_Tabas(Tf,dt):
    df = pd.read_excel (r'TABAS.xlsx')
    baseacc=np.array(df)
    baseacc = np.insert(baseacc, 0, np.array([0,0.009]), axis=0)
    NT=floor(Tf/dt)
    t=np.linspace(0,Tf,NT)
    new_acc = np.interp(t, baseacc[:,0], baseacc[:,1])*9.80665
    return new_acc

def ddxg_Loma(Tf,dt):
    df = pd.read_excel (r'Loma.xlsx')
    baseacc=np.array(df)
    baseacc = np.insert(baseacc, 0, np.array([0,0]), axis=0)
    NT=floor(Tf/dt)
    t=np.linspace(0,Tf,NT)
    new_acc = np.interp(t, baseacc[:,0], baseacc[:,1])*9.80665
    return new_acc

def Newmark_coefficients(dt):
    alpha=0.25
    delta=0.5

    a0=1/(alpha*(dt**2))
    a1=delta/(alpha*dt)
    a2=1/(alpha*dt)
    a3=(1/(2*alpha))-1
    a4=(delta/alpha)-1
    a5=(dt/2)*((delta/alpha)-2)
    a6=dt*(1-delta)
    a7=delta*dt

    return a0,a1,a2,a3,a4,a5,a6,a7

################################################################################

def Newmark_initialize(ndof,a0,a1,M,C,K,NT,DI,VI):
    
    nlength=len(atleast_1d(M)) 

    KH=zeros((ndof,ndof),float)

    KH=K+a0*M+a1*C

    if(nlength==1):

        U=zeros(NT,float)
        Ud=zeros(NT,float)
        Udd=zeros(NT,float)

        U[0]=DI
        Ud[0]=VI
            
    else:

        U=zeros((ndof,NT),float)
        Ud=zeros((ndof,NT),float)
        Udd=zeros((ndof,NT),float)

        U[:,0]=DI
        Ud[:,0]=VI
    return U,Ud,Udd,KH

def Newmark_PID(M,alpha,K,r,kp,ki,kd,VI,DI,dt,Tf,ndof,ddxg):
    """
    input

    M = mass matrix
    alpha = hyperparameter (for zeta configuration in order to obtain damping matrix)
    K = stiffness matrix
    


    VI = initial velocity
    DI = initial displacement

      dt = time step
      NT = number of time points
    ndof = number of degrees of freedom

    output

      U = displacement
     Ud = velocity
    Udd = acceleration

    """
#uncontrolled sys    
    a0,a1,a2,a3,a4,a5,a6,a7=Newmark_coefficients(dt)
    NT=floor(Tf/dt)
    C=dampparams(alpha,M,K)
    U,Ud,Udd,KH=Newmark_initialize(ndof,a0,a1,M,C,K,NT,DI,VI)
    t=np.linspace(0,Tf,NT)
    f1=0
    f2=0
    
    for i in range (1,NT):
        R=np.array([f1-M[0,0]*ddxg[i-1],f2-M[1,1]*ddxg[i-1]])
        V1=(a1*U[:,i-1]+a4*Ud[:,i-1]+a5*Udd[:,i-1])
        V2=(a0*U[:,i-1]+a2*Ud[:,i-1]+a3*Udd[:,i-1])
        CV=dot(C,V1)
        MA=dot(M,V2)
        FH=R+MA+CV

#  solve for displacements

        if(ndof>1):
            Un = linalg.solve(KH, FH)
        else:
            Un=FH/KH

        Uddn=a0*(Un-U[:,i-1])-a2*Ud[:,i-1]-a3*Udd[:,i-1]
        Udn=Ud[:,i-1]+a6*Udd[:,i-1]+a7*Uddn


        U[:,i]=Un
        Ud[:,i]=Udn
        Udd[:,i]=Uddn
        
#### control system
    e=0
    f2c=0
    f1c=0
    integral_prev=0
    eprev = 0
    error=[]
    f1cont=[]
    error.append(0)
    f1cont.append(0)
    DDxg=[]
    DDxg.append(ddxg)
    deltast1c=[]
    deltast2c=[]
    deltast1c.append((f1c+f2c-3*M[1,1]*ddxg)/K[1,1])
    deltast2c.append((f1c+2*f2c-4*M[1,1]*ddxg)/K[1,1]);
    Uc,Udc,Uddc,KHc=Newmark_initialize(ndof,a0,a1,M,C,K,NT,DI,VI)
    for i in range (1,NT):
        e = r-Uc[1,i-1]; #sensor is located on second mass
        integeral = integral_prev + e*dt
        f1c = kp*e + ki*integeral + kd*(error[i-2]-4*eprev+3*e)/(2*dt)
        Rc=np.array([f1c-M[0,0]*ddxg[i-1],f2c-M[1,1]*ddxg[i-1]])
        V1c=(a1*Uc[:,i-1]+a4*Udc[:,i-1]+a5*Uddc[:,i-1])
        V2c=(a0*Uc[:,i-1]+a2*Udc[:,i-1]+a3*Uddc[:,i-1])
        CVc=dot(C,V1c)
        MAc=dot(M,V2c)
        FHc=Rc+MAc+CVc

#  solve for displacements

        if(ndof>1):
            Unc = linalg.solve(KHc, FHc)
        else:
            Unc=FHc/KHc

        Uddnc=a0*(Unc-Uc[:,i-1])-a2*Udc[:,i-1]-a3*Uddc[:,i-1]
        Udnc=Udc[:,i-1]+a6*Uddc[:,i-1]+a7*Uddnc

        eprev = e;
        integral_prev = integeral
        Uc[:,i]=Unc
        Udc[:,i]=Udnc
        Uddc[:,i]=Uddnc
        error.append(e)
        f1cont.append(f1c)
        deltast1c.append((f1c+f2c-3*M[1,1]*ddxg)/K[1,1])
        deltast2c.append((f1c+2*f2c-4*M[1,1]*ddxg)/K[1,1])
        DDxg.append(ddxg)
        
#  delta static
    deltast1 = (f1+f2-3*M[1,1]*ddxg)/K[1,1]
    deltast2 = (f1+2*f2-4*M[1,1]*ddxg)/K[1,1]

    mpl.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 16
    plt.rcParams['axes.linewidth'] = 2
    #Displacements
    fig,axs = plt.subplots(2, 1)
    fig.set_size_inches(18.5, 10.5, forward=True)
    axs[0].step(t,U[0,:], c='gray',linewidth=2.5)
    axs[0].step(t,Uc[0,:],c='black',linewidth=2.5)
    axs[0].axhline(y = deltast1[0], color = 'black', linestyle = '-.')
    axs[0].axhline(y = r, color = 'black', linestyle = '--')
   # axs[0].set_title('Displacement (m)',fontsize=16)
    axs[0].set_xlim(0,Tf)
    axs[0].legend(['Uncontrolled','PID',r'Static equilibrium $x_{1st}$'])#'x1 = x1* - deltast1','x1c = x1c* - deltast1c'
    axs[0].set_ylabel(r'$ x_1$ (m)',fontsize=18)#\/\/
    axs[0].label_outer()
    axs[1].step(t,U[1,:],c='gray',linewidth=2.5)
    axs[1].step(t,Uc[1,:],c='black',linewidth=2.5)
    axs[1].axhline(y = deltast2[0], color = 'black', linestyle = '-.')
    axs[1].axhline(y = r, color = 'black', linestyle = '--')
    axs[1].set_xlim(0,Tf)
    axs[1].set_xlabel('time (sec)',fontsize=18)
    axs[1].set_ylabel(r'$ x_2$ (m)',fontsize=18)
    axs[1].legend(['Uncontrolled','PID',r'Static equilibrium $x_{2st}$'])#'x2 = x2* - deltast2','x2c = x2c* - deltast2c'
    fig.savefig('Displacemet.svg', transparent=False, bbox_inches='tight')
    
    #Velocities
    fig1,axs1 = plt.subplots(2, 1)
    fig1.set_size_inches(18.5, 10.5, forward=True)
    axs1[0].step(t,Ud[0,:])
    axs1[0].step(t,Udc[0,:])
    axs1[0].set_xlim(0,Tf)
    axs1[0].axhline(y = 0, color = 'black', linestyle = '--')
    axs1[0].set_title('Velocity (m/s)',fontsize=16)
    axs1[0].set_ylabel("x'1 (m/s)",fontsize=14)
    axs1[0].legend(["x'1","x'1c"])
    axs1[0].label_outer()
    axs1[1].step(t,Ud[1,:])
    axs1[1].step(t,Udc[1,:])
    axs1[1].axhline(y = 0, color = 'black', linestyle = '--')
    axs1[1].set_xlim(0,Tf)
    axs1[1].set_ylabel("x'2 (m/s)",fontsize=14)
    axs1[1].set_xlabel('time (sec)',fontsize=14)
    axs1[1].legend(["x'2","x'2c"])
    
    #Accelerations
    fig2,axs2 = plt.subplots(2, 1)
    fig2.set_size_inches(18.5, 10.5, forward=True)
    axs2[0].step(t,Udd[0,:])
    axs2[0].step(t,Uddc[0,:])
    axs2[0].set_xlim(0,Tf)
    axs2[0].axhline(y = 0, color = 'black', linestyle = '--')
    axs2[0].set_title('Acceleration (m/s^2)',fontsize=16)
    axs2[0].set_ylabel('x"1 (m/s^2)',fontsize=14)
    axs2[0].legend(['x"1','x"1c'])
    axs2[0].label_outer()
    axs2[1].step(t,Udd[1,:])
    axs2[1].step(t,Uddc[1,:])
    axs2[1].set_xlim(0,Tf)
    axs2[1].axhline(y = 0, color = 'black', linestyle = '--')
    axs2[1].set_ylabel('x"2 (m/s^2)',fontsize=14)
    axs2[1].set_xlabel('time (sec)',fontsize=14)
    axs2[1].legend(['x"2','x"2c'])
    
    #Error and control force
    fig3,axs3 = plt.subplots(2, 1)
    fig3.set_size_inches(18.5, 10.5, forward=True)
    axs3[0].step(t,error,c='black',linewidth=2.5)
    axs3[0].set_xlim(0,Tf)
    axs3[0].axhline(y = 0, color = 'black', linestyle = '--')
   # axs3[0].set_title('Error Signal',fontsize=16)
    axs3[0].set_ylabel('Error e(t)' ,fontsize=18)
    axs3[0].label_outer()
    axs3[1].step(t,f1cont,c='black',linewidth=2.5)
    axs3[1].set_xlim(0,Tf)
    axs3[1].axhline(y = f1cont[NT-1], color = 'black', linestyle = '--')
    axs3[1].set_ylabel('Control force (N)',fontsize=18)
    axs3[1].set_xlabel('time (sec)',fontsize=18)
    fig3.savefig('controlforce.svg', transparent=False, bbox_inches='tight')
    
    #disturbances
    fig4,axs4 = plt.subplots(1, 1)
    fig4.set_size_inches(18.5, 5.5, forward=True)
    axs4.step(t,ddxg)
    axs4.set_xlim(0,Tf)
    axs4.set_title('Base Acceleration',fontsize=16)
    axs4.set_ylabel('ddxg (m/s^2)',fontsize=14)
    axs4.label_outer()
    plt.show()

