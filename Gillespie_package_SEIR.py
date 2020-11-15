import numpy as np

class SEIR_package_model_Gillespie:
    def __init__(self, Ndays, occupancy, incidence_frac, infection_probs, wp_contact_rate, \
                 td, alphaT, alphaC, sigma, gamma, cross_contact_rate, parcels_per_day, del_hlife, virus_hlife):
        self.Ndays = Ndays  #scalar
        if len(occupancy) > 1:
            self.occupancy = np.array(occupancy)
        else:
            self.occupancy = np.ones(self.Ndays)*occupancy[0]
        if len(incidence_frac) > 1:
            self.incidence = np.array(incidence_frac)
        else:
            self.incidence = np.ones(self.Ndays)*incidence_frac[0]
        if len(parcels_per_day) > 1:
            self.parcels = np.array(parcels_per_day)
        else:
            self.parcels = np.ones(self.Ndays)*parcels_per_day[0]
        self.infection_prob = infection_probs #always a scalar
        if len(wp_contact_rate) > 1:
            self.wp_contact_rate = np.array(wp_contact_rate)
        else:
            self.wp_contact_rate = np.ones(self.Ndays)*wp_contact_rate[0]
        self.td = td
        self.alphaT = alphaT
        self.alphaC = alphaC
        self.sigma = sigma
        self.gamma = gamma
        self.phi = cross_contact_rate
        self.del_rate = np.log(2)/del_hlife
        self.ster_rate = np.log(2)/virus_hlife
        self.cmat = np.zeros((3,3))
        self.cmat[0,0] = self.td
        self.cmat[0,1] = self.phi*self.td
        self.cmat[0,2] = self.phi*self.td
        self.cmat[1,0] = self.cmat[0,1]
        self.cmat[1,1] = 1
        self.cmat[1,2] = self.phi
        self.cmat[2,0] = self.cmat[0,2]
        self.cmat[2,1] = self.cmat[1,2]
        self.cmat[2,2] = 1
    
    def get_lambda_parts(self,N0,s,e,i,ip,ipaft,time):
        p = np.zeros(15)
        it = int(time)
        bmat = -self.wp_contact_rate[it]*np.log(1-self.infection_prob)*self.occupancy[it]*self.cmat
        ivect = np.dot(bmat,i)/np.sum(N0)
        for k in np.arange(3):  
            p[k] = (self.incidence[it] + ivect[k])*s[k] #infection events
            p[k+3] = self.sigma*e[k]    #prodromal events
            p[k+6] = self.gamma*i[k]    #recovery events
        p[0] += self.alphaC*ip*s[0]/N0[0]     #contracted from package
        p[9] = self.alphaT*(self.parcels[it]-ip-ipaft)*i[1]/N0[1]     #transmission to package
        p[10] = self.alphaT*(self.parcels[it]-ip-ipaft)*i[0]/N0[0]     #transmission to package
        p[11] = self.del_rate*ip 
        p[12] = self.del_rate*ipaft               #package removed
        p[13] = self.ster_rate*ip               #package sterilised
        p[14] = self.ster_rate*ipaft
        
        return p
    
    def generate_dt(self, N0,s,e,i,ip,ipaft,time):   
        u = np.random.rand(1)
        dt = np.ceil(time) - time
        cdf = 0
        dtime = 0
        while cdf < -np.log(1-u) and time + dtime < self.Ndays:
            lambda0_p = self.get_lambda_parts(N0,s,e,i,ip,ipaft,time + dtime)
            lambda0 = np.sum(lambda0_p)
            cdf += lambda0*dt
            dtime += dt
            dt = 1
        
        if cdf >= -np.log(1-u):
            dtime -= (cdf + np.log(1-u))/lambda0
            return dtime
        else:  #no events remaining
            return -1

    def generate_event(self, N0,s,e,i,ip,ipaft,time):
        u = np.random.rand(1)
        lam_p = self.get_lambda_parts(N0,s,e,i,ip,ipaft,time)
        cdf = np.cumsum(lam_p)/np.sum(lam_p)
        event = np.sum(u > cdf)
        
        return event
    
    def update_state(self,s,e,i,ip,ipaft,dp,rp,event):
        sh = s.copy()
        eh = e.copy()
        ih = i.copy()
        if event < 3:
            sh[event] -= 1
            eh[event] += 1
        else:
            if event < 6:
                eh[event-3] -= 1
                ih[event-3] += 1
            else:
                if event < 9:
                    ih[event-6] -= 1
                else:
                    if event == 9:
                        ip += 1
                    else:
                        if event == 10:
                            ipaft += 1
                        else:
                            if event == 11 or event==13:
                                ip -= 1
                            else:
                                ipaft -= 1
                            if event == 11 or event==12:
                                dp += 1
                            else:
                                rp += 1
                
        return sh,eh,ih,ip,ipaft,dp,rp
    
    def run_sim(self,N0,S0,E0,I0,ip,ipaft):
        s = S0.copy()
        e = E0.copy()
        i = I0.copy()
        dp = 0
        rp = 0
        state_space = np.zeros((13,1))
        state_space[0:3,0] = s
        state_space[3:6,0] = e
        state_space[6:9,0] = i
        state_space[9,0] = ip
        state_space[10,0] = ipaft
        state_space[11,0] = dp
        state_space[12,0] = rp
        time = np.zeros(1)
        k = 1
        t = 0
        while t < self.Ndays:
            dt = self.generate_dt(N0,s,e,i,ip,ipaft,t)
            if dt >= 0:  #steady state
                t += dt
                event = self.generate_event(N0,s,e,i,ip,ipaft,t)
                s,e,i,ip,ipaft,dp,rp = self.update_state(s,e,i,ip,ipaft,dp,rp,event)
            else:
                t = self.Ndays
            #print(dt,event)
            state_space = np.append(state_space,np.zeros((13,1)),axis=1)
            state_space[0:3,k] = s
            state_space[3:6,k] = e
            state_space[6:9,k] = i
            state_space[9,k] = ip
            state_space[10,k] = ipaft
            state_space[11,k] = dp
            state_space[12,k] = rp
            time = np.append(time,t)
            k += 1
            
        return time, state_space[0:3,:], state_space[3:6,:], state_space[6:9,:], state_space[9:13,:]



    def fill_SEI_disease_free_Jacobian(self,N0,i_time=0):
        J = np.zeros((7,7))
        bmat = -self.wp_contact_rate[i_time]*np.log(1-self.infection_prob)*self.occupancy[i_time]*self.cmat
        Np = self.parcels[i_time]
        for i in np.arange(3):
            J[i,i] -= self.sigma
            J[3+i,i] += self.sigma
            J[3+i,3+i] -= self.gamma
            for j in np.arange(3):
                J[i,3+j] += bmat[i,j]*N0[i]/np.sum(N0)
        J[6,4] += self.alphaT*Np*np.sum(N0)/N0[1]
        J[0,6] += self.alphaC/N0[0]
        J[6,6] -= (self.del_rate + self.ster_rate)

        return J

    
    def run_deterministic_sim(self,N0,ip,dt):
        Nt = int(self.Ndays/dt)
        S = np.zeros((3,Nt))
        E = np.zeros((3,Nt))
        I = np.zeros((3,Nt))
        R = np.zeros((3,Nt))
        Parcels = np.zeros((5,Nt))
        
        F = np.array(N0)/np.sum(N0)
        #need to work out step equations, as it is non-linear
        S[:,0] = (np.array(N0) - np.array(ip))/np.sum(N0)
        E[:,0] = np.array(ip)/np.sum(N0)
        Parcels[0,0] = self.parcels[0]
        
        Xold = np.append(np.append(np.append(S[:,0],E[:,0]),I[:,0]),Parcels[:3,0])
        for n in np.arange(1,Nt):
            day = int(n*dt)
            Xprev = np.append(np.append(np.append(S[:,n-1],E[:,n-1]),I[:,n-1]),Parcels[:3,n-1])
            Xnew = Xold.copy()
            f = self.model_function(Xprev,Xnew,N0,day,dt)
            count = 0
            while np.linalg.norm(f)/dt > 1E-10:
                J = self.fill_Jacobian(Xprev,Xnew,N0,day,dt)
                Xnew = Xold - np.linalg.solve(J,f)
                f = self.model_function(Xprev,Xnew,N0,day,dt)
                Xold = Xnew.copy()
                count += 1
            S[:,n] = Xnew[:3]
            E[:,n] = Xnew[3:6]
            I[:,n] = Xnew[6:9]
            Parcels[:3,n] = Xnew[9:12]
            R[:,n] = F - S[:,n] - E[:,n] - I[:,n]
            Parcels[3,n] = Parcels[3,n-1] + dt*(self.del_rate)*(Parcels[1,n] + Parcels[2,n])
            Parcels[4,n] = Parcels[4,n-1] + dt*(self.ster_rate)*(Parcels[1,n] + Parcels[2,n])
    
        return S,E,I,R,Parcels

    def model_function(self,Xprev,Xnew,N0,i_time,dt):
        f = Xnew - Xprev
        bmat = -self.wp_contact_rate[i_time]*np.log(1-self.infection_prob)*self.occupancy[i_time]*self.cmat
        f[:3] += dt*Xnew[:3]*(np.dot(bmat,Xnew[6:9]) + self.incidence[i_time])
        f[0] += dt*self.alphaC*Xnew[0]*(Xnew[10]/N0[0])
        f[3] += -dt*self.alphaC*Xnew[0]*(Xnew[10]/N0[0])
        f[3:6] += -dt*Xnew[:3]*(np.dot(bmat,Xnew[6:9]) + self.incidence[i_time]) + dt*self.sigma*Xnew[3:6] 
        f[6:9] += -dt*self.sigma*Xnew[3:6] + dt*self.gamma*Xnew[6:9]
        f[9] += dt*Xnew[9]*(self.alphaT*np.sum(N0)*(Xnew[7]/N0[1] + Xnew[6]/N0[0]) + self.del_rate) - self.parcels[i_time]*dt
        f[10] += -dt*Xnew[9]*(np.sum(N0)*self.alphaT*Xnew[7]/N0[1]) + dt*(self.del_rate + self.ster_rate)*Xnew[10]
        f[11] += -dt*Xnew[9]*(np.sum(N0)*self.alphaT*Xnew[6]/N0[0]) + dt*(self.del_rate + self.ster_rate)*Xnew[11]

        return f

    def fill_Jacobian(self,Xprev,Xnew,N0,i_time,dt):
        J = np.zeros((12,12))
        bmat = -self.wp_contact_rate[i_time]*np.log(1-self.infection_prob)*self.occupancy[i_time]*self.cmat
        for i in np.arange(12):
            J[i,i] += 1
        J[0,0] += dt*self.alphaC*Xnew[10]/N0[0]
        J[0,10] += dt*self.alphaC*Xnew[0]/N0[0]
        J[3,0] += -dt*self.alphaC*Xnew[10]/N0[0]
        J[3,10] += -dt*self.alphaC*Xnew[0]/N0[0]
        J[9,6] += dt*self.alphaT*np.sum(N0)*Xnew[9]/N0[0]
        J[9,7] += dt*self.alphaT*np.sum(N0)*Xnew[9]/N0[1]
        J[9,9] += dt*(self.alphaT*np.sum(N0)*(Xnew[7]/N0[1] + Xnew[6]/N0[0]) + self.del_rate)
        J[10,7] += -dt*self.alphaT*np.sum(N0)*Xnew[9]/N0[1]
        J[10,9] += -dt*self.alphaT*np.sum(N0)*Xnew[7]/N0[1]
        J[10,10] += dt*(self.del_rate + self.ster_rate)
        J[11,6] += -dt*self.alphaT*np.sum(N0)*Xnew[9]/N0[0]
        J[11,9] += -dt*self.alphaT*np.sum(N0)*Xnew[6]/N0[0]
        J[11,11] += dt*(self.del_rate + self.ster_rate)
        for i in np.arange(3):
            J[i,i] += dt*(np.dot(bmat,Xnew[6:9])[i] + self.incidence[i_time])
            J[3+i,i] += -dt*(np.dot(bmat,Xnew[6:9])[i] + self.incidence[i_time])
            J[3+i,3+i] += dt*self.sigma
            J[6+i,3+i] += -dt*self.sigma
            J[6+i,6+i] += dt*self.gamma
            for j in np.arange(3):
                J[i,6+j] += dt*Xnew[i]*bmat[i,j]
                J[3+i,6+j] += -dt*Xnew[i]*bmat[i,j]

        return J