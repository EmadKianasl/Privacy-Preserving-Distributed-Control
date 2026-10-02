############SECURE CONTROLLED SYS #########################entire system##########################################
#########PID controller + Uncontrolled system

from numpy import zeros,dot,atleast_1d
import numpy as np
import pandas as pd
from scipy import linalg
from matplotlib import pyplot as plt
from math import sqrt,floor
import time
import numpy as np
import shamirproc
import shamirscheme as sss
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from tabulate import tabulate
from decimal import Decimal
import matplotlib as mpl
from scipy.io import savemat
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

def Newmark_PID(F,server,t,scale,M,alpha,K,r,kp,ki,kd,VI,DI,dt,Tf,ndof,ddxg):

#uncontrolled sys    
    a0,a1,a2,a3,a4,a5,a6,a7=Newmark_coefficients(dt)
    NT=floor(Tf/dt)
    C=dampparams(alpha,M,K)
    U,Ud,Udd,KH=Newmark_initialize(ndof,a0,a1,M,C,K,NT,DI,VI)
    taim=np.linspace(0,Tf,NT)
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
        
        
#### insecure control system
    e=0
    f2c=0
    f1c=0
    integral_prev=0
    eprev = 0
    error=[]
    f1cont=[]
    error.append(0)
    f1cont.append(0)
  #  DDxg=[]
   # DDxg.append(ddxg)
   # deltast1c=[]
    #deltast2c=[]
   # deltast1c.append((f1c+f2c-3*M[1,1]*ddxg)/K[1,1])
   # deltast2c.append((f1c+2*f2c-4*M[1,1]*ddxg)/K[1,1]);
    Uc,Udc,Uddc,KHc=Newmark_initialize(ndof,a0,a1,M,C,K,NT,DI,VI)
    tic=time.time()
    for i in range (1,NT):
        e = r-Uc[1,i-1]; #sensor is located on the second mass
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
      #  deltast1c.append((f1c+f2c-3*M[1,1]*ddxg)/K[1,1])
      #  deltast2c.append((f1c+2*f2c-4*M[1,1]*ddxg)/K[1,1])
    toc=time.time()    
    
    print('field cardinality : {0}'.format(F.cardinality()))
#### secure control system
    es=sss.share(F,0,t,server,scale)
    f2cs=0
    f1cs=0
    integral_prevs=sss.share(F,0,t,server,scale)
    eprevs = sss.share(F,0,t,server,scale)
    errors=[]
    f1conts=[]
    f1conts_shares=[]
    errorsh = []
    errorsh.append(eprevs)
    errors.append(0)
    f1conts.append(0)
    #DDxgs=[]
    #DDxgs.append(ddxg)
    #deltast1cs=[]
    #deltast2cs=[]
    #deltast1cs.append((f1cs+f2cs-3*M[1,1]*ddxg)/K[1,1])
    #deltast2cs.append((f1cs+2*f2cs-4*M[1,1]*ddxg)/K[1,1])
    Ucs,Udcs,Uddcs,KHcs=Newmark_initialize(ndof,a0,a1,M,C,K,NT,DI,VI)
    rs = sss.share(F,r,t,server,scale)
    Ucs2 = sss.share(F,-Ucs[1,0],t,server,scale)
    rec=sss.basispoly(F,server)
    dts=sss.share(F,dt,t,server,scale)
    invdts=sss.share(F,1/(2*dt),t,server,scale)
    kp=sss.share(F,kp,t,server,scale)
    ki=sss.share(F,ki,t,server,scale)
    kd=sss.share(F,kd,t,server,scale)
    tics=time.time()
    for i in range (1,NT):
        es = shamirproc.add(F,rs,Ucs2,t,server,scale) #sensor is located on the second mass
        integerals = shamirproc.add(F,integral_prevs, shamirproc.mult(F,es,dts,t,server,scale),t,server,scale)
        f11cs = shamirproc.mult(F,kp,es,t,server,scale) #proportional
        f12cs = shamirproc.mult(F,ki,integerals,t,server,scale) #Integral
        f13cs = shamirproc.mult(F,shamirproc.mult(F,kd,shamirproc.add(F,3*es,
                                                                      shamirproc.add(F,-4*eprevs,errorsh[i-2],t,server,scale),t,server,scale)
                                                  ,t,server,scale),invdts,t,server,scale) #derivative
        f1cs = shamirproc.add(F,f11cs,shamirproc.add(F,f12cs,f13cs,t,server,scale),t,server,scale)  
        f1conts_shares.append(f1cs)
        
        Rcs=np.array([sss.rec(F,f1cs,rec,scale)-M[0,0]*ddxg[i-1],f2cs-M[1,1]*ddxg[i-1]])
        V1cs=(a1*Ucs[:,i-1]+a4*Udcs[:,i-1]+a5*Uddcs[:,i-1])
        V2cs=(a0*Ucs[:,i-1]+a2*Udcs[:,i-1]+a3*Uddcs[:,i-1])
        CVcs=dot(C,V1cs)
        MAcs=dot(M,V2cs)
        FHcs=Rcs+MAcs+CVcs
    
#  solve for displacements

        if(ndof>1):
            Uncs = linalg.solve(KHcs, FHcs)
        else:
            Uncs=FHcs/KHcs

        Uddncs=a0*(Uncs-Ucs[:,i-1])-a2*Udcs[:,i-1]-a3*Uddcs[:,i-1]
        Udncs=Udcs[:,i-1]+a6*Uddcs[:,i-1]+a7*Uddncs
        Ucs2=sss.share(F,-Uncs[1],t,server,scale)
        
        eprevs=es
        integral_prevs = integerals
        Ucs[:,i]=Uncs
        Udcs[:,i]=Udncs
        Uddcs[:,i]=Uddncs
        errorsh.append(es)
        errors.append(sss.rec(F,es,rec,scale))
        f1conts.append(sss.rec(F,f1cs,rec,scale))
        #deltast1cs.append((f1cs+f2cs-3*M[1,1]*ddxg)/K[1,1])
       # deltast2cs.append((f1cs+2*f2cs-4*M[1,1]*ddxg)/K[1,1])
    tocs=time.time()
    print('running time of insecure controller: {0} s'.format(toc-tic))
    print('running time of secure controller: {0} s'.format(tocs-tics))
    print('running time of secure /running time of insecure : {0}'.format((tocs-tics)/(toc-tic)))
    
    disp1=U[0,:]
    disp2=U[1,:]
    vel1=Ud[0,:]
    vel2=Ud[1,:]
    acc1=Udd[0,:]
    acc2=Udd[1,:]
    
    disp1c=Uc[0,:]
    disp2c=Uc[1,:]
    vel1c=Udc[0,:]
    vel2c=Udc[1,:]
    acc1c=Uddc[0,:]
    acc2c=Uddc[1,:]
    
    disp1cs=Ucs[0,:]
    disp2cs=Ucs[1,:]
    vel1cs=Udcs[0,:]
    vel2cs=Udcs[1,:]
    acc1cs=Uddcs[0,:]
    acc2cs=Uddcs[1,:]
    
    #Frequency Responses
    freq= np.fft.fftfreq(NT,dt)
    freq= np.fft.fftshift(freq)
    disp1_fft = np.fft.fft(disp1)*dt
    disp1_fft= np.fft.fftshift(disp1_fft)
    disp2_fft = np.fft.fft(disp2)*dt
    disp2_fft= np.fft.fftshift(disp2_fft)
    vel1_fft = np.fft.fft(vel1)*dt
    vel1_fft= np.fft.fftshift(vel1_fft)
    vel2_fft = np.fft.fft(vel2)*dt
    vel2_fft= np.fft.fftshift(vel2_fft)
    acc1_fft = np.fft.fft(acc1)*dt
    acc1_fft= np.fft.fftshift(acc1_fft)
    acc2_fft = np.fft.fft(acc2)*dt 
    acc2_fft= np.fft.fftshift(acc2_fft)
    disp1c_fft = np.fft.fft(disp1c)*dt
    disp1c_fft= np.fft.fftshift(disp1c_fft)
    disp2c_fft = np.fft.fft(disp2c)*dt
    disp2c_fft= np.fft.fftshift(disp2c_fft)
    vel1c_fft = np.fft.fft(vel1c)*dt
    vel1c_fft= np.fft.fftshift(vel1c_fft)
    vel2c_fft = np.fft.fft(vel2c)*dt
    vel2c_fft= np.fft.fftshift(vel2c_fft)
    acc1c_fft = np.fft.fft(acc1c)*dt
    acc1c_fft= np.fft.fftshift(acc1c_fft)
    acc2c_fft = np.fft.fft(acc2c)*dt 
    acc2c_fft= np.fft.fftshift(acc2c_fft)
    disp1cs_fft = np.fft.fft(disp1cs)*dt
    disp1cs_fft= np.fft.fftshift(disp1cs_fft)
    disp2cs_fft = np.fft.fft(disp2cs)*dt
    disp2cs_fft= np.fft.fftshift(disp2cs_fft)
    vel1cs_fft = np.fft.fft(vel1cs)*dt
    vel1cs_fft= np.fft.fftshift(vel1cs_fft)
    vel2cs_fft = np.fft.fft(vel2cs)*dt
    vel2cs_fft= np.fft.fftshift(vel2cs_fft)
    acc1cs_fft = np.fft.fft(acc1cs)*dt
    acc1cs_fft= np.fft.fftshift(acc1cs_fft)
    acc2cs_fft = np.fft.fft(acc2cs)*dt   
    acc2cs_fft= np.fft.fftshift(acc2cs_fft)
    error_fft = np.fft.fft(error)*dt
    error_fft= np.fft.fftshift(error_fft)
    errors_fft = np.fft.fft(errors)*dt
    errors_fft= np.fft.fftshift(errors_fft)
    f1cont_fft = np.fft.fft(f1cont)*dt
    f1cont_fft= np.fft.fftshift(f1cont_fft)
    f1conts_fft = np.fft.fft(f1conts)*dt
    f1conts_fft= np.fft.fftshift(f1conts_fft)
    
    #passing params
    f1conts_shares=np.array(f1conts_shares,dtype='uint64')
    Newmark_PID.f1conts_shares=f1conts_shares
    Newmark_PID.errors=errors
    Newmark_PID.f1conts=f1conts
    Newmark_PID.error=error
    Newmark_PID.f1cont=f1cont
    Newmark_PID.disp1=disp1
    Newmark_PID.disp2=disp2
    Newmark_PID.vel1=vel1
    Newmark_PID.vel2=vel2
    Newmark_PID.acc1=acc1
    Newmark_PID.acc2=acc2
    Newmark_PID.disp1c=disp1c
    Newmark_PID.disp2c=disp2c
    Newmark_PID.vel1c=vel1c
    Newmark_PID.vel2c=vel2c
    Newmark_PID.acc1c=acc1c
    Newmark_PID.acc2c=acc2c
    Newmark_PID.disp1cs=disp1cs
    Newmark_PID.disp2cs=disp2cs
    Newmark_PID.vel1cs=vel1cs
    Newmark_PID.vel2cs=vel2cs
    Newmark_PID.acc1cs=acc1cs
    Newmark_PID.acc2cs=acc2cs
    Newmark_PID.taim=taim
    taim2=np.delete(taim, -1)
    Newmark_PID.taim2=taim2
    
    Newmark_PID.freq= freq
    Newmark_PID.disp1_fft =disp1_fft 
    Newmark_PID.disp2_fft =disp2_fft
    Newmark_PID.vel1_fft =vel1_fft 
    Newmark_PID.vel2_fft = vel2_fft 
    Newmark_PID.acc1_fft =acc1_fft 
    Newmark_PID.acc2_fft =acc2_fft
    Newmark_PID.disp1c_fft = disp1c_fft 
    Newmark_PID.disp2c_fft =disp2c_fft
    Newmark_PID.vel1c_fft =vel1c_fft
    Newmark_PID.vel2c_fft =vel2c_fft
    Newmark_PID.acc1c_fft =acc1c_fft
    Newmark_PID.acc2c_fft =acc2c_fft
    Newmark_PID.disp1cs_fft =disp1cs_fft
    Newmark_PID.disp2cs_fft=disp2cs_fft
    Newmark_PID. vel1cs_fft = vel1cs_fft
    Newmark_PID.vel2cs_fft = vel2cs_fft
    Newmark_PID.acc1cs_fft = acc1cs_fft
    Newmark_PID.acc2cs_fft =acc2cs_fft  
    Newmark_PID.error_fft =error_fft
    Newmark_PID.errors_fft =errors_fft
    Newmark_PID.f1cont_fft = f1cont_fft 
    Newmark_PID.f1conts_fft =f1conts_fft 
    
    
    
    
    results={
    'f1conts_shares':f1conts_shares,
    'errors':errors,
    'error':error,
    'f1conts':f1conts ,
    'f1cont':f1cont,
    'disp1':disp1,
    'disp2':disp2,
    'vel1':vel1,
    'vel2':vel2,
    'acc1':acc1,
    'acc2':acc2,
    'disp1c':disp1c,
    'disp2c':disp2c,
    'vel1c':vel1c,
    'vel2c':vel2c,
    'acc1c':acc1c,
    'acc2c':acc2c,
    'disp1cs':disp1cs,
    'disp2cs':disp2cs,
    'vel1cs':vel1cs,
    'vel2cs':vel2cs,
    'acc1cs':acc1cs,
    'acc2cs':acc2cs,
    'taim':taim,
    'taim2':taim2,
    'freq':freq,
    'disp1_fft':np.abs(disp1_fft) ,
    'disp2_fft' :np.abs(disp2_fft),
    'vel1_fft' :np.abs(vel1_fft),
    'vel2_fft' : np.abs(vel2_fft) ,
    'acc1_fft' :np.abs(acc1_fft) ,
    'acc2_fft' :np.abs(acc2_fft),
    'disp1c_fft' : np.abs(disp1c_fft) ,
    'disp2c_fft' :np.abs(disp2c_fft),
    'vel1c_fft' :np.abs(vel1c_fft),
    'vel2c_fft' :np.abs(vel2c_fft),
    'acc1c_fft' :np.abs(acc1c_fft),
    'acc2c_fft' :np.abs(acc2c_fft),
    'disp1cs_fft' :np.abs(disp1cs_fft),
    'disp2cs_fft':np.abs(disp2cs_fft),
    'vel1cs_fft' : np.abs(vel1cs_fft),
    'vel2cs_fft' : np.abs(vel2cs_fft),
    'acc1cs_fft' : np.abs(acc1cs_fft),
    'acc2cs_fft' :np.abs(acc2cs_fft) , 
    'error_fft' :np.abs(error_fft),
    'errors_fft' :np.abs(errors_fft),
    'f1cont_fft' :np.abs( f1cont_fft) ,
    'f1conts_fft' :np.abs(f1conts_fft)  
    }
    savemat("results.mat",results)
    #f1conts_shares=f1conts_shares.tolist()
    #savemat("f1conts_shares.mat", f1conts_shares)
    #savemat("errors.mat", errors)
    #savemat("error.mat", error)
    #savemat("f1conts.mat", f1conts)
    #savemat("f1cont.mat", f1cont)
    #savemat("disp1.mat", disp1)
    #savemat("disp2.mat", disp2)
    #savemat("vel1.mat", vel1)
    #savemat("vel2.mat", vel2)
    #savemat("acc1.mat", acc1)
    #savemat("acc2.mat", acc2)
    #savemat("disp1c.mat", disp1c)
    #savemat("disp2c.mat", disp2c)
    #savemat("vel1c.mat", vel1c)
    #savemat("vel2c.mat", vel2c)
    #savemat("acc1c.mat", acc1c)
    #savemat("acc2c.mat", acc2c)
    #savemat("disp1cs.mat", disp1cs)
    #savemat("disp2cs.mat", disp2cs)
    #savemat("vel1cs.mat", vel1cs)
    #savemat("vel2cs.mat", vel2cs)
   # savemat("acc1cs.mat", acc1cs)
  #  savemat("acc2cs.mat", acc2cs)
  #  savemat("taim.mat", taim)
    
    #signal energy and power Uncontrolled
    x1E=0
    for i in abs(np.array(U[0,:])):
        x1E+=i**2
    x1P='%.3E' % Decimal(str(x1E/NT))
    
    x2E=0
    for i in abs(np.array(U[1,:])):
        x2E+=i**2
    x2P='%.3E' % Decimal(str(x2E/NT))
    
    dx1E=0
    for i in abs(np.array(Ud[0,:])):
        dx1E+=i**2
    dx1P='%.3E' % Decimal(str(dx1E/NT))
    
    dx2E=0
    for i in abs(np.array(Ud[1,:])):
        dx2E+=i**2
    dx2P='%.3E' % Decimal(str(dx2E/NT))    
    
    ddx1E=0
    for i in abs(np.array(Udd[0,:])):
        ddx1E+=i**2
    ddx1P='%.3E' % Decimal(str(ddx1E/NT))

    ddx2E=0
    for i in abs(np.array(Udd[1,:])):
        ddx2E+=i**2
    ddx2P='%.3E' % Decimal(str(ddx2E/NT)) 
    
    ddxgE=0
    for i in abs(np.array(ddxg)):
        ddxgE+=i**2
    ddxgP='%.3E' % Decimal(str(ddxgE/NT)) 
    
    #signal energy and power controlled-insecure
    
    x1cE=0
    for i in abs(np.array(Uc[0,:])):
        x1cE+=i**2
    x1cP='%.3E' % Decimal(str(x1cE/NT))
    
    x2cE=0
    for i in abs(np.array(Uc[1,:])):
        x2cE+=i**2
    x2cP='%.3E' % Decimal(str(x2cE/NT))
    
    dx1cE=0
    for i in abs(np.array(Udc[0,:])):
        dx1cE+=i**2
    dx1cP='%.3E' % Decimal(str(dx1cE/NT))
    
    dx2cE=0
    for i in abs(np.array(Udc[1,:])):
        dx2cE+=i**2
    dx2cP='%.3E' % Decimal(str(dx2cE/NT))    
    
    ddx1cE=0
    for i in abs(np.array(Uddc[0,:])):
        ddx1cE+=i**2
    ddx1cP='%.3E' % Decimal(str(ddx1cE/NT))

    ddx2cE=0
    for i in abs(np.array(Uddc[1,:])):
        ddx2cE+=i**2
    ddx2cP='%.3E' % Decimal(str(ddx2cE/NT))
    
    ecE=0
    for i in abs(np.array(error)):
        ecE+=i**2
    ecP='%.3E' % Decimal(str(ecE/NT))
    
    fcE=0
    for i in abs(np.array(f1cont)):
        fcE+=i**2
    fcP='%.3E' % Decimal(str(fcE/NT))
    
    
    #signal energy and power controlled-secure
    
    x1csE=0
    for i in abs(np.array(Ucs[0,:])):
        x1csE+=i**2
    x1csP='%.3E' % Decimal(str(x1csE/NT))
    
    x2csE=0
    for i in abs(np.array(Ucs[1,:])):
        x2csE+=i**2
    x2csP='%.3E' % Decimal(str(x2csE/NT))
    
    dx1csE=0
    for i in abs(np.array(Udcs[0,:])):
        dx1csE+=i**2
    dx1csP='%.3E' % Decimal(str(dx1csE/NT))
    
    dx2csE=0
    for i in abs(np.array(Udcs[1,:])):
        dx2csE+=i**2
    dx2csP='%.3E' % Decimal(str(dx2csE/NT))    
    
    ddx1csE=0
    for i in abs(np.array(Uddcs[0,:])):
        ddx1csE+=i**2
    ddx1csP='%.3E' % Decimal(str(ddx1csE/NT))

    ddx2csE=0
    for i in abs(np.array(Uddcs[1,:])):
        ddx2csE+=i**2
    ddx2csP='%.3E' % Decimal(str(ddx2csE/NT))
   
    ecsE=0
    for i in abs(np.array(errors)):
        ecsE+=i**2
    ecsP='%.3E' % Decimal(str(ecsE/NT))
    
    fcsE=0
    for i in abs(np.array(f1conts)):
        fcsE+=i**2
    fcsP='%.3E' % Decimal(str(fcsE/NT))
 
    
    ##differences
    x1diffrence = abs(np.array(Uc[0,:])-np.array(Ucs[0,:]))
    dx1diffrence = abs(np.array(Udc[0,:])-np.array(Udcs[0,:]))
    ddx1diffrence = abs(np.array(Uddc[0,:])-np.array(Uddcs[0,:]))
    x2diffrence = abs(np.array(Uc[1,:])-np.array(Ucs[1,:]))
    dx2diffrence = abs(np.array(Udc[1,:])-np.array(Udcs[1,:]))
    ddx2diffrence = abs(np.array(Uddc[1,:])-np.array(Uddcs[1,:]))
    errordiffrence = abs(np.array(error)-np.array(errors))
    forcediffrence = abs(np.array(f1cont)-np.array(f1conts))
    
    #table of defferences
    sx1=0
    for i in x1diffrence:
        sx1+=i**2
    IAEx1= '%.3E' % Decimal(str(sum(x1diffrence)/NT))
    ISEx1= '%.3E' % Decimal(str(sx1/NT))

    sx2=0
    for i in x2diffrence:
        sx2+=i**2
    IAEx2= '%.3E' % Decimal(str(sum(x2diffrence)/NT))
    ISEx2= '%.3E' % Decimal(str(sx2/NT))
    
    sdx1=0
    for i in dx1diffrence:
        sdx1+=i**2
    IAEdx1= '%.3E' % Decimal(str(sum(dx1diffrence)/NT))
    ISEdx1= '%.3E' % Decimal(str(sdx1/NT))

    sdx2=0
    for i in dx2diffrence:
        sdx2+=i**2
    IAEdx2= '%.3E' % Decimal(str(sum(dx2diffrence)/NT))
    ISEdx2= '%.3E' % Decimal(str(sdx2/NT))
    
    sddx1=0
    for i in ddx1diffrence:
        sddx1+=i**2
    IAEddx1= '%.3E' % Decimal(str(sum(ddx1diffrence)/NT))
    ISEddx1= '%.3E' % Decimal(str(sddx1/NT))

    sddx2=0
    for i in ddx2diffrence:
        sddx2+=i**2
    IAEddx2= '%.3E' % Decimal(str(sum(ddx2diffrence)/NT))
    ISEddx2= '%.3E' % Decimal(str(sddx2/NT))
    
    se=0
    for i in errordiffrence:
        se+=i**2
    IAEse= '%.3E' % Decimal(str(sum(errordiffrence)/NT))
    ISEse= '%.3E' % Decimal(str(se/NT))
    
    sf=0
    for i in forcediffrence:
        sf+=i**2
    IAEsf= '%.3E' % Decimal(str(sum(forcediffrence)/NT))
    ISEsf= '%.3E' % Decimal(str(sf/NT))
    
    table = [["x1",IAEx1,ISEx1,x1P,x1cP,x1csP]
             ,["x2",IAEx2,ISEx2,x2P,x2cP,x2csP]
             ,["dx1/dt",IAEdx1,ISEdx1,dx1P,dx1cP,dx1csP]
            ,["dx2/dt",IAEdx2,ISEdx2,dx2P,dx2cP,dx2csP]
            ,['ddx1/dt',IAEddx1,ISEddx1,ddx1P,ddx1cP,ddx1csP]
            ,['ddx2/dt',IAEddx2,ISEddx2,ddx2P,ddx2cP,ddx2csP]
            ,["PID error signal",IAEse,ISEse,0,ecP,ecsP]
            ,["PID force",IAEsf,ISEsf,0,fcP,fcsP]
            ,['x"g',0,0,ddxgP,0,0]]

    headers = ['Signal',"Integral Absolute Error(IAE)\nsum[abs(insecureX-secureX)]/samples", "Integral Squared Error(ISE)\nsum[abs(sqrt(insecureX-secureX))]/samples",'Open Loop','Closed Loop/PID','Closed Loop/MPC-PID']

    print(tabulate(table, headers, tablefmt="pretty",disable_numparse=True))
    fold=open('table.txt','w')
    fold.write(tabulate(table))
    fold.close()
##  percentage error    
   # pex1=(x1diffrence/np.array(Uc[0,:]))*100
    #pex2=(x2diffrence/np.array(Uc[1,:]))*100
    #pedx1=(dx1diffrence/np.array(Udc[0,:]))*100
    #pedx2=(dx2diffrence/np.array(Udc[1,:]))*100
    #peddx1=abs(ddx1diffrence/np.array(Uddc[0,:]))*100
    #peddx2=abs(ddx2diffrence/np.array(Uddc[1,:]))*100
    #peerror=abs(errordiffrence/np.array(error))*100
   # peforce=abs(forcediffrence/np.array(f1cont))*100
#  delta static
    #deltast1 = (f1+f2-3*M[1,1]*ddxg)/K[1,1]
    #deltast2 = (f1+2*f2-4*M[1,1]*ddxg)/K[1,1]

    mpl.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 16
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams["mathtext.fontset"]='dejavuserif'
    space = 3
    dash_len = 3
    #Displacements
    fig,axs = plt.subplots(2,1)
    fig.set_size_inches(18.5, 10.5, forward=True)
  #  fig.suptitle('Displacement (m)', fontsize=16)
    
    axs[0].step(taim,U[0,:], color='gray')
    axs[0].step(taim,Uc[0,:], color='orange')
    axs[0].step(taim,Ucs[0,:], color='darkblue',linestyle = '--',dashes=(dash_len, space))
   # axs[0,0].axhline(y = deltast1, color = 'brown', linestyle = '--')
    #axs[0,0].axhline(y = r, color = 'black', linestyle = '--')
    axs[0].set_xlim(0,Tf)
    axs[0].legend(['Uncontrolled','PID','SMPC-PID'])
    axs[0].set_xlabel('Time (sec)')
    axs[0].set_ylabel(r'Displacement $x_1$ (m)')
   # axs[0,0].set_title('Insecure')
    axs[0].grid()
    axs[0].label_outer()
    
    #axs[0,1].step(taim,U[0,:], color='red')
    #axs[0,1].step(taim,Ucs[0,:], color='blue')
  #  axs[0,1].axhline(y = deltast1, color = 'brown', linestyle = '--')
   # axs[0,1].axhline(y = r, color = 'black', linestyle = '--')
    #axs[0,1].set_title('Secure')
    #axs[0,1].set_xlim(0,Tf)
    #axs[0,1].grid()
    #axs[0,1].legend(['Uncontrolled','Secure Distributed PID Controller'])

    
   # axs[0,2].step(taim,x1diffrence, color='purple')
   # axs[0,2].grid()
   # axs[0,2].set_title('insecure and secure overall difference')
   # axs[0,2].set_xlim(0,Tf)
    
    axs[1].step(taim,U[1,:], color='gray')
    axs[1].step(taim,Uc[1,:], color='orange')
    axs[1].step(taim,Ucs[1,:], color='darkblue',linestyle ='--',dashes=(dash_len, space))
   # axs[1,0].axhline(y = deltast2, color = 'brown', linestyle = '--')
    #axs[1,0].axhline(y = r, color = 'black', linestyle = '--')
    axs[1].set_xlim(0,Tf)
    axs[1].set_xlabel('Time (sec)')
    axs[1].set_ylabel(r'Displacement $x_2$ (m)')
    axs[1].grid()
    axs[1].legend(['Uncontrolled','PID','SMPC-PID'])
    fig.savefig('Displacements2.svg', transparent=False, bbox_inches='tight')
   # 
    #axs[1,1].step(taim,U[1,:], color='red')
    #axs[1,1].step(taim,Ucs[1,:], color='blue')
   # axs[1,1].axhline(y = deltast2, color = 'brown', linestyle = '--')
    #axs[1,1].axhline(y = r, color = 'black', linestyle = '--')
    #axs[1,1].set_xlim(0,Tf)
    #axs[1,1].set_xlabel('time (sec)',fontsize=14)
    #axs[1,1].grid()
    #axs[1,1].legend(['x2','x2c_secure'])

  #  axs[1,2].step(taim,x2diffrence, color='purple')
   # axs[1,2].grid()
    #axs[1,2].set_xlim(0,Tf)

    #Velocities
    fig1,axs1 = plt.subplots(2, 1)
    fig1.set_size_inches(18.5, 10.5, forward=True)
   # fig1.suptitle('Velocity (m/s)', fontsize=16)
    
    axs1[0].step(taim,Ud[0,:], color='gray')
    axs1[0].step(taim,Udc[0,:], color='orange')
    axs1[0].step(taim,Udcs[0,:],color='darkblue',linestyle ='--',dashes=(dash_len, space))
    #axs1[0].set_title('Insecure')
    axs1[0].set_xlim(0,Tf)
    axs1[0].legend(['Uncontrolled','PID','SMPC-PID'])
    axs1[0].set_xlabel('Time (sec)')
    axs1[0].set_ylabel(r"Velocity $\.x_1$ (m/s)")
   # axs1[0].axhline(y = 0, color = 'black', linestyle = '--')
    axs1[0].grid()
    axs1[0].label_outer()
    
  #  axs1[0,1].step(taim,Ud[0,:], color='red')
   # axs1[0,1].step(taim,Udcs[0,:], color='blue')
    #axs1[0,1].set_title('Secure')
    #axs1[0,1].set_xlim(0,Tf)
    #axs1[0,1].legend(["x'1","x'1c_secure"])
    #axs1[0,1].grid()
   # axs1[0,1].axhline(y = 0, color = 'black', linestyle = '--')

    
  #  axs1[0,2].step(taim,dx1diffrence, color='purple')
   # axs1[0,2].axhline(y = 0, color = 'black', linestyle = '--')
   # axs1[0,2].set_title('insecure and secure overall difference')
   # axs1[0,2].grid()
   # axs1[0,2].set_xlim(0,Tf)
    
    axs1[1].step(taim,Ud[1,:], color='gray')
    axs1[1].step(taim,Udc[1,:], color='orange')
    axs1[1].step(taim,Udcs[1,:],color='darkblue',linestyle ='--',dashes=(dash_len, space))
    axs1[1].set_xlim(0,Tf)
    axs1[1].set_xlabel('Time (sec)')
    axs1[1].set_ylabel(r"Velocity $\.x_2$ (m/s)")
   # axs1[1].axhline(y = 0, color = 'black', linestyle = '--')
    axs1[1].grid()
    axs1[1].legend(['Uncontrolled','PID','SMPC-PID'])
    fig1.savefig('Velocities_2.svg', transparent=False, bbox_inches='tight')

   # axs1[1,1].step(taim,Ud[1,:], color='red')
   # axs1[1,1].step(taim,Udcs[1,:], color='blue')
    #axs1[1,1].set_xlim(0,Tf)
    #axs1[1,1].legend(["x'2","x'2c_secure"])
  #  axs1[1,1].axhline(y = 0, color = 'black', linestyle = '--')
    #axs1[1,1].grid()
    #axs1[1,1].set_xlabel('time (sec)',fontsize=14)
    
    #axs1[1,2].step(taim,dx2diffrence, color='purple')
    #axs1[1,2].axhline(y = 0, color = 'black', linestyle = '--')
    #axs1[1,2].grid()
    #axs1[1,2].set_xlim(0,Tf)
    
    
    #Accelerations
    fig2,axs2 = plt.subplots(2,1)
    fig2.set_size_inches(18.5, 10.5, forward=True)    
  #  fig2.suptitle('َAcceleration (m/s^2)', fontsize=16)
    
    axs2[0].step(taim,Udd[0,:], color='gray')
    axs2[0].step(taim,Uddc[0,:], color='orange')
    axs2[0].step(taim,Uddcs[0,:],color='darkblue',linestyle ='--',dashes=(dash_len, space))
   # axs2[0,0].axhline(y = 0, color = 'black', linestyle = '--')
   # axs2[0].set_title('Insecure')
    axs2[0].set_xlim(0,Tf)
    axs2[0].legend(['Uncontrolled','PID','SMPC-PID'])
    axs2[0].set_xlabel('Time (sec)')
    axs2[0].set_ylabel(r"Acceleration $\ddot x_1\ (m/s^2)$")
    axs2[0].grid()
    axs2[0].label_outer()
    
    
  #  axs2[0,1].step(taim,Udd[0,:], color='red')
   # axs2[0,1].step(taim,Uddcs[0,:], color='blue')
   # axs2[0,1].axhline(y = 0, color = 'black', linestyle = '--')
    #axs2[0,1].set_title('Secure')
    #axs2[0,1].set_xlim(0,Tf)
    #axs2[0,1].grid()
    #axs2[0,1].legend(['x"1','x"1c_secure'])
    
  #  axs2[0,2].step(taim,ddx1diffrence, color='purple')
  #  axs2[0,2].axhline(y = 0, color = 'black', linestyle = '--')
  #  axs2[0,2].set_title('insecure and secure overall difference')
   # axs2[0,2].grid()
    #axs2[0,2].set_xlim(0,Tf)
    
    axs2[1].step(taim,Udd[1,:], color='gray')
    axs2[1].step(taim,Uddc[1,:], color='orange')
    axs2[1].step(taim,Uddcs[1,:],color='darkblue',linestyle ='--',dashes=(dash_len, space))
   # axs2[1].axhline(y = 0, color = 'black', linestyle = '--')
    axs2[1].set_xlim(0,Tf)
    axs2[1].set_xlabel('Time (sec)')
    axs2[1].set_ylabel(r"Acceleration $\ddot x_2\ (m/s^2)$")
    axs2[1].grid()
    axs2[1].legend(['Uncontrolled','PID','SMPC-PID'])
    fig2.savefig('Accelerations2.svg', transparent=False, bbox_inches='tight')
    
    #axs2[1,1].step(taim,Udd[1,:], color='red')
    #axs2[1,1].step(taim,Uddcs[1,:], color='blue')
   ## axs2[1,1].axhline(y = 0, color = 'black', linestyle = '--')
    #axs2[1,1].set_xlim(0,Tf)
    #axs2[1,1].set_xlabel('time (sec)',fontsize=14)
    #axs2[1,1].grid()
    #axs2[1,1].legend(['x"2','x"2c_secure'])    
    
    #axs2[1,2].step(taim,ddx2diffrence, color='purple')
    #axs2[1,2].axhline(y = 0, color = 'black', linestyle = '--')
    #axs2[1,2].grid()
    #axs2[1,2].set_xlim(0,Tf)

    #Error and control force
    fig3,axs3 = plt.subplots(2, 1)
    fig3.set_size_inches(18.5, 10.5, forward=True)
   # fig3.suptitle('error signal and control force', fontsize=16)
    axs3[0].step(taim,error, color='orange')
    axs3[0].step(taim,errors,color='darkblue',linestyle ='--',dashes=(dash_len, space))
    #axs3[0].set_title('Insecure')
    axs3[0].set_xlim(0,Tf)
    axs3[0].set_ylabel('Error (m)')
    axs3[0].legend(['PID','SMPC-PID '])  
   # axs3[0].axhline(y = 0, color = 'black', linestyle = '--')
    axs3[0].grid()
    axs3[0].label_outer()
    
  #  axs3[0,1].step(taim,errors, color='blue')
   # axs3[0,1].set_title('Secure')
    #axs3[0,1].set_xlim(0,Tf)
    #axs3[0,1].legend(['error signal_secure'])
    #axs3[0,1].grid()
  #  axs3[0,1].axhline(y = 0, color = 'black', linestyle = '--')
    
    #axs3[0,2].step(taim,errordiffrence, color='purple')
    #axs3[0,2].axhline(y = 0, color = 'black', linestyle = '--')
    #axs3[0,2].set_title('insecure and secure overall difference ')
    #axs3[0,2].grid()
    #axs3[0,2].set_xlim(0,Tf)
    
    axs3[1].step(taim,f1cont, color='orange')
    axs3[1].step(taim,f1conts,color='darkblue',linestyle ='--',dashes=(dash_len, space))
    axs3[1].set_xlim(0,Tf)
    axs3[1].set_ylabel('Control force (N)')
    axs3[1].set_xlabel('Time (sec)')
   # axs3[1].axhline(y = f1conts[-1], color = 'black', linestyle = '--')
    axs3[1].grid()
    axs3[1].legend(['PID','SMPC-PID '])
    fig3.savefig('Errors and forces 2.svg', transparent=False, bbox_inches='tight')
    
    #axs3[1,1].step(taim,f1conts, color='blue')
    #axs3[1,1].set_xlim(0,Tf)
    #axs3[1,1].set_xlabel('time (sec)',fontsize=14)
    #axs3[1,1].axhline(y = f1conts[-1], color = 'black', linestyle = '--')
    #axs3[1,1].grid()
    #axs3[1,1].legend(['control force_secure'])
    
    #axs3[1,2].step(taim,forcediffrence, color='purple')
    #axs3[1,2].axhline(y = 0, color = 'black', linestyle = '--')
    #axs3[1,2].grid()
    #axs3[1,2].set_xlim(0,Tf)
    

    #disturbances
    fig4,axs4 = plt.subplots(1, 1)
    fig4.set_size_inches(18.5, 5.5, forward=True)
    axs4.step(taim,ddxg, color='blue')
    axs4.set_xlim(0,Tf)
    axs4.set_title('Base Acceleration',fontsize=16)
    axs4.set_ylabel('ddxg (m/s^2)',fontsize=14)
    axs4.grid()
    axs4.label_outer()
    

    #server actual computation
    fig5, axs5 = plt.subplots()
    fig5.set_size_inches(9.5, 5, forward=True)
    for i in range(server):
        axs5.step(taim2,f1conts_shares[:,i])
              
    axs5.set_xlim(0,5)
    axs5.set_ylabel('Distributed control force (N)')
    axs5.set_xlabel('Time (sec)')
    leg=[]
    for i in range(server):
        leg.append('Server {0}'.format(str(i+1)))
    axs5.legend(leg, loc='upper right')
    fig5.savefig('shared force 2.svg', transparent=False, bbox_inches='tight')
    
    ###############Frequency Responses plots###################
    
    #displacements
    fig6, axs6 = plt.subplots(2,1)
    fig6.set_size_inches(18.5, 10.5, forward=True)
   # fig1.suptitle('Velocity (m/s)', fontsize=16)
    
    axs6[0].semilogx(freq,np.abs(disp1_fft), color='gray')
    axs6[0].semilogx(freq,np.abs(disp1c_fft), color='orange')
    axs6[0].semilogx(freq,np.abs(disp1cs_fft),color='darkblue',linestyle ='--',dashes=(dash_len, space))
    #axs1[0].set_title('Insecure')
    axs6[0].set_xlim(0,10)
    axs6[0].legend(['Uncontrolled','PID','SMPC-PID'])
    axs6[0].set_xlabel('Frequency (Hz)')
    axs6[0].set_ylabel("Spectral Displacement")
   # axs1[0].axhline(y = 0, color = 'black', linestyle = '--')
    axs6[0].grid()
    axs6[0].label_outer()
    
    
    axs6[1].semilogx(freq,np.abs(disp2_fft), color='gray')
    axs6[1].semilogx(freq,np.abs(disp2c_fft), color='orange')
    axs6[1].semilogx(freq,np.abs(disp2cs_fft),color='darkblue',linestyle ='--',dashes=(dash_len, space))
    axs6[1].set_xlim(0,10)
    axs6[1].set_xlabel('Frequency (Hz)')
    axs6[1].set_ylabel("Spectral Displacement")
   # axs1[1].axhline(y = 0, color = 'black', linestyle = '--')
    axs6[1].grid()
    axs6[1].legend(['Uncontrolled','PID','SMPC-PID'])

    
    
    ##velocities
    fig7,axs7 = plt.subplots(2, 1)
    fig7.set_size_inches(18.5, 10.5, forward=True)
   # fig1.suptitle('Velocity (m/s)', fontsize=16)
    
    axs7[0].semilogx(freq,np.abs(vel1_fft), color='gray')
    axs7[0].semilogx(freq,np.abs(vel1c_fft), color='orange')
    axs7[0].semilogx(freq,np.abs(vel1cs_fft),color='darkblue',linestyle ='--',dashes=(dash_len, space))
    #axs1[0].set_title('Insecure')
    axs7[0].set_xlim(0,10)
    axs7[0].legend(['Uncontrolled','PID','SMPC-PID'])
    axs7[0].set_xlabel('Frequency (Hz)')
    axs7[0].set_ylabel("Spectral Velocity")
   # axs1[0].axhline(y = 0, color = 'black', linestyle = '--')
    axs7[0].grid()
    axs7[0].label_outer()
        
    axs7[1].semilogx(freq,np.abs(vel2_fft), color='gray')
    axs7[1].semilogx(freq, np.abs(vel2c_fft), color='orange')
    axs7[1].semilogx(freq, np.abs(vel2cs_fft),color='darkblue',linestyle ='--',dashes=(dash_len, space))
    axs7[1].set_xlim(0,10)
    axs7[1].set_xlabel('Frequency (Hz)')
    axs7[1].set_ylabel("Spectral Velocity")
   # axs1[1].axhline(y = 0, color = 'black', linestyle = '--')
    axs7[1].grid()
    axs1[1].legend(['Uncontrolled','PID','SMPC-PID'])
    
     #Accelerations
    fig8,axs8 = plt.subplots(2,1)
    fig8.set_size_inches(18.5, 10.5, forward=True)    
  #  fig2.suptitle('َAcceleration (m/s^2)', fontsize=16)
    
    axs8[0].semilogx(freq,np.abs(acc1_fft), color='gray')
    axs8[0].semilogx(freq,np.abs(acc1c_fft), color='orange')
    axs8[0].semilogx(freq,np.abs(acc1cs_fft),color='darkblue',linestyle ='--',dashes=(dash_len, space))
   # axs2[0,0].axhline(y = 0, color = 'black', linestyle = '--')
   # axs2[0].set_title('Insecure')
    axs8[0].set_xlim(0,10)
    axs8[0].legend(['Uncontrolled','PID','SMPC-PID'])
    axs8[0].set_xlabel('Frequency (Hz)')
    axs8[0].set_ylabel("Spectral Acceleration")
    axs8[0].grid()
    axs8[0].label_outer()
    
    axs8[1].semilogx(freq,np.abs(acc2_fft), color='gray')
    axs8[1].semilogx(freq,np.abs(acc2c_fft), color='orange')
    axs8[1].semilogx(freq,np.abs(acc2cs_fft),color='darkblue',linestyle ='--',dashes=(dash_len, space))
   # axs2[1].axhline(y = 0, color = 'black', linestyle = '--')
    axs8[1].set_xlim(0,10)
    axs8[1].set_xlabel('Frequency (Hz)')
    axs8[1].set_ylabel("Spectral Acceleration")
    axs8[1].grid()
    axs8[1].legend(['Uncontrolled','PID','SMPC-PID'])
  #  fig8.savefig('Accelerations2.svg', transparent=False, bbox_inches='tight')
    
   #Error and control force
    fig9,axs9 = plt.subplots(2, 1)
    fig9.set_size_inches(18.5, 10.5, forward=True)
   # fig3.suptitle('error signal and control force', fontsize=16)
    axs9[0].semilogx(freq,np.abs(error_fft), color='orange')
    axs9[0].semilogx(freq,np.abs(errors_fft),color='darkblue',linestyle ='--',dashes=(dash_len, space))
    #axs3[0].set_title('Insecure')
    axs9[0].set_xlim(0,10)
    axs9[0].set_ylabel('Spectral Error')
    axs9[0].legend(['PID','SMPC-PID '])  
   # axs3[0].axhline(y = 0, color = 'black', linestyle = '--')
    axs9[0].grid()
    axs9[0].label_outer()
    
    axs9[1].semilogx(freq,np.abs(f1cont_fft), color='orange')
    axs9[1].semilogx(freq,np.abs(f1conts_fft),color='darkblue',linestyle ='--',dashes=(dash_len, space))
    axs9[1].set_xlim(0,10)
    axs9[1].set_ylabel('Spectral Control force')
    axs9[1].set_xlabel('Frequency (Hz)')
   # axs3[1].axhline(y = f1conts[-1], color = 'black', linestyle = '--')
    axs9[1].grid()
    axs9[1].legend(['PID','SMPC-PID '])
  #  fig3.savefig('Errors and forces 2.svg', transparent=False, bbox_inches='tight')   