import numpy as np

def SIRsim(beta,td,sigma,gamma,FOI,ND,NL,NO,alpha1,alpha2,pert,dt=0.1,Nt=1000):
    """Run SIR sim for parameters: 
    beta = baseline transmission rate for workplace (per day)
    td = fraction of time drivers spend at depot
    sigma = inverse of latent period
    gamma = revocery rate (inverse of infectious period)
    FOI = external force of infection
    ND = number of drivers
    NL = number of loaders
    NO = number of other/office staff
    alpha1 = rate of fomite transmission on to packages
    alpha2 = rate of fomite transmission from infectious packages
    dt = timestep
    Nt = number of timesteps
    Simulation uses a backwards Euler formulation with a Newton-Raphson solve for each timestep
    """
    S = np.zeros((4,Nt+1))
    E = np.zeros((3,Nt+1))
    I = np.zeros((4,Nt+1))
    R = np.zeros((4,Nt+1))
    N = np.array([ND,NL,NO])/(NL+NO+ND)
    #need to work out step equations, as it is non-linear
    S[0:3,0] = N
    S[1,0] -= pert
    S[3,0] = 1
    E[1,0] = pert
    
    Xold = np.append(np.append(S[:,0],E[:,0]),I[:,0])
    for n in np.arange(1,Nt+1):
        Xprev = np.append(np.append(S[:,n-1],E[:,n-1]),I[:,n-1])
        Xnew = Xold.copy()
        f = model_function(Xprev,Xnew,beta,td,sigma,gamma,FOI,N,alpha1,alpha2,dt,Nt)
        while np.linalg.norm(f)/dt > 1E-12:
            J = fill_Jacobian(Xprev,Xnew,beta,td,sigma,gamma,FOI,N,alpha1,alpha2,dt,Nt)
            Xnew = Xold - np.linalg.solve(J,f)
            f = model_function(Xprev,Xnew,beta,td,sigma,gamma,FOI,N,alpha1,alpha2,dt,Nt)
            Xold = Xnew.copy()
        S[:,n] = Xnew[:4]
        E[:,n] = Xnew[4:7]
        I[:,n] = Xnew[7:]
        R[0:3,n] = N - S[0:3,n] - E[0:3,n] - I[0:3,n]
        R[3,n] = R[3,n-1] + dt*I[3,n] + dt*alpha2*S[3,n]*I[0,n]/S[0,0]
    
    return S,E,I,R

def model_function(Xprev,Xnew,nu,td,sigma,gamma,FOI,N,alpha1,alpha2,dt,Nt):
    t = np.array([td,1,1])
    beta = nu*np.outer(t,t)/(N[0]*td + N[1] + N[2])
    f = Xnew - Xprev
    f[:3] += dt*Xnew[:3]*(np.dot(beta,Xnew[7:10]) + FOI)
    f[0] += dt*alpha2*Xnew[0]*(Xnew[10]/N[0])
    f[3] += dt*Xnew[3]*(alpha1*Xnew[8]/N[1] + 1) - dt
    f[4:7] += -dt*Xnew[:3]*(np.dot(beta,Xnew[7:10]) + FOI) + dt*sigma*Xnew[4:7] 
    f[4] += - dt*Xnew[0]*alpha2*(Xnew[10]/N[0])
    f[7:10] += -dt*sigma*Xnew[4:7] + dt*gamma*Xnew[7:10]
    f[10] += -dt*Xnew[3]*(alpha1*Xnew[8]/N[1]) + dt*Xnew[10]
    
    return f

def fill_Jacobian(Xprev,Xnew,nu,td,sigma,gamma,FOI,N,alpha1,alpha2,dt,Nt):
    J = np.zeros((11,11))
    t = np.array([td,1,1])
    beta = nu*np.outer(t,t)/(N[0]*td + N[1] + N[2])
    for i in np.arange(11):
        J[i,i] = 1
    J[0,0] += dt*alpha2*Xnew[10]/N[0]
    J[0,10] += dt*alpha2*Xnew[0]/N[0]
    J[4,0] -= dt*alpha2*Xnew[10]/N[0]
    J[4,10] -= dt*alpha2*Xnew[0]/N[0]
    for i in np.arange(3):
        J[i,i] += dt*(np.dot(beta[i,:],Xnew[7:10]) + FOI)
        J[4+i,i] += -dt*(np.dot(beta[i,:],Xnew[7:10]) + FOI)
        J[4+i,4+i] += dt*sigma
        J[7+i,4+i] += -dt*sigma
        J[7+i,7+i] += dt*gamma
        for j in np.arange(3):
            J[i,7+j] += dt*Xnew[i]*beta[i,j]
            J[4+i,7+j] += -dt*Xnew[i]*beta[i,j]
    J[3,3] += dt*(alpha1*Xnew[8]/N[1] + 1)
    J[3,8] += dt*alpha1*Xnew[3]/N[1]
    J[10,3] += -dt*alpha1*Xnew[8]/N[1]
    J[10,8] += -dt*alpha1*Xnew[3]/N[1]
    J[10,10] += dt
    
    return J