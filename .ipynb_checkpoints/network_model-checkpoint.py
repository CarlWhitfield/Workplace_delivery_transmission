import numpy as np

PAsymp = 0.3
Pisolate = 0.5
IsolationLength = 14

def get_gamma_params(mean, std):
    k = mean**2/std**2
    theta = std**2/mean
    return k,theta

#these params need to be set
k_l, theta_l = get_gamma_params(5, 2)
k_p, theta_p = get_gamma_params(2, 1)
k_s, theta_s = get_gamma_params(7, 5)
k_a, theta_a = get_gamma_params(5, 2)

def generate_latent_period():
    time = np.random.gamma(k_l, theta_l)
    day = int(np.round(time))
    if day == 0:
        day = 1
    return day
    
def generate_prodromal_period():
    time = np.random.gamma(k_p, theta_p)
    day = int(np.round(time))
    if day == 0:
        day = 1
    return day
    
def generate_symptomatic_period(Asymp):
    if Asymp:
        time = np.random.gamma(k_a, theta_a)
    else:
        time = np.random.gamma(k_s, theta_s)
    day = int(np.round(time))
    if day == 0:
        day = 1
    return day
    
class Node:
    def __init__(self, node_type = 'd'):
        self.type = node_type
        self.state = 0      #s = 0, e = 1, p = 2, i = 3, r = 4
        self.isolating = 0
    
    def increment_state(self):
        if self.state < 4:
            self.state += 1
        else:
            self.state = 0
    
    def update_status(self, day):
        if self.state > 0:
            dayr = day - self.infection_day
            if dayr == self.latent_period:
                self.increment_state()
            if dayr - self.latent_period == self.prodromal_period:
                self.increment_state()
                if self.asymp == 0:
                    self.isolating = np.random.binomial(1,Pisolate)
            if dayr - self.latent_period - self.prodromal_period == self.symptomatic_period:
                self.increment_state()
            if dayr - self.latent_period - self.prodromal_period == IsolationLength:
                self.isolating = 0
                
    def infect_node(self, infection_day, infected_by):
        if self.state == 0:
            self.increment_state()
            self.latent_period = generate_latent_period()
            self.prodromal_period = generate_prodromal_period()
            self.asymp = np.random.binomial(1,PAsymp)
            self.infectivity = self.asymp*0.5 + (1-self.asymp) 
            self.symptomatic_period = generate_symptomatic_period(self.asymp)
            self.infected_by = infected_by
            self.infection_day = infection_day
        
class Workplace_Network:
    def __init__(self, N, wd_occupancy, sat_occupancy, onsite_shifts_per_day = 1,\
                 driver_pairs=True, fixed_pairings = False):
        #N is vector of numbers in each role: 0 = drivers, 1 = loaders, 2 = office
        self.ND = int(N[0])
        self.NL = int(N[1])
        self.NO = int(N[2])
        self.nodes = []
        for n in np.arange(self.ND):
            self.nodes.append(Node('d'))
        for n in np.arange(self.NL):
            self.nodes.append(Node('l'))
        for n in np.arange(self.NO):
            self.nodes.append(Node('o'))  
        self.infectious_nodes = []
        self.pairs = driver_pairs
        self.fixed_pairs = fixed_pairings
        self.day_rate = np.array(5*[wd_occupancy] + [sat_occupancy] + [0])
        self.os_shifts = onsite_shifts_per_day
    
    
    def run_simulation(self, Ndays, contacts_pd, x_contact_factor, inf_prob, 
                       pair_inf_prob, starting_infections=1, EFOI = 0):
        self.beta = inf_prob
        self.kappa = contacts_pd
        self.phi = x_contact_factor
        self.alpha = pair_inf_prob
        self.states = np.zeros((Ndays,5,3))
        self.Nisolating = np.zeros((Ndays,3))
        SI = np.random.choice(len(self.nodes),starting_infections)
        for i in SI:
            self.nodes[i].infect_node(0,-1)
        
        for n in np.arange(Ndays):
            new_infs = self.run_single_day(n,EFOI)
            self.states[n,:,:], self.Nisolating[n,:] = self.count_states()
    
    def run_single_day(self, day_no, external_infection_rate):
        new_infections = []
        self.infectious_nodes = []
        drivers_not_isolating = []
        loaders_not_isolating = []
        office_not_isolating = []
        #update status of nodes
        for i in np.arange(len(self.nodes)):
            self.nodes[i].update_status(day_no)
            #check which nodes are at work
            if self.nodes[i].isolating == 0:
                if self.nodes[i].type == 'd':
                    drivers_not_isolating.append(i)
                if self.nodes[i].type == 'l':
                    loaders_not_isolating.append(i)
                if self.nodes[i].type == 'o':
                    office_not_isolating.append(i)    
                #check which nodes are infectious (p and i both infectious states)
                if self.nodes[i].state == 2 or self.nodes[i].state == 3:
                    self.infectious_nodes.append(i)  
            
            i += 1
        
        #do new external introductions
        Ni = np.random.binomial(self.ND+self.NL+self.NO, external_infection_rate)
        nodes_ext_infected = np.random.choice(self.ND+self.NL+self.NO,Ni,replace = False)
        for i in nodes_ext_infected:
            if self.nodes[i].state ==0:
                new_infections.append(self.nodes[i])
                self.nodes[i].infect_node(day_no,-1)  #only works if node is susceptible, -1 indicates external infection
            
        if self.day_rate[day_no%7] == 0 or len(self.infectious_nodes)==0:
            return new_infections
        
        #draw employees who are at work
        if self.pairs:
            d_pairs = self.draw_driver_pairs(drivers_not_isolating, self.day_rate[day_no%7])
            drivers = np.array(d_pairs).flatten()
        else:
            NDdraw = int(np.round(self.day_rate[day_no%7]*self.ND))
            if NDdraw < len(drivers_not_isolating):
                drivers = np.random.choice(drivers_not_isolating, NDdraw, replace = False)
            else:
                drivers = np.array(drivers_not_isolating)
        NLdraw = int(np.round(self.day_rate[day_no%7]*self.NL))
        if NLdraw < len(loaders_not_isolating):
            loaders = np.random.choice(loaders_not_isolating, NLdraw, replace = False)
        else:
            loaders = np.array(loaders_not_isolating)
        NOdraw = int(np.round(self.day_rate[day_no%7]*self.NO))
        if NOdraw < len(office_not_isolating):
            office = np.random.choice(office_not_isolating, NOdraw, replace = False)
        else:
            office = np.array(office_not_isolating)
        
        NDh = len(drivers)
        NLh = len(loaders)
        NOh = len(office)
        #find out which employees are infectious and draw infectious contacts
        #infect only if susceptible
        #split office and loaders into shifts randomly
        for i in self.infectious_nodes:
            #generate contacts
            Ncons = np.random.poisson(self.kappa)
            #infectious contacts
            Ninf_cons = np.random.binomial(Ncons, self.beta)
            #draw contacts at random
            if self.nodes[i].type == 'd':
                NDh -= 1
                Probs = NDh*[1] + (NLh + NOh)*[self.phi]
            if self.nodes[i].type == 'l':
                NLh -= 1
                Probs = NDh*[self.phi] + NLh*[1] + NOh*[self.phi]
            if self.nodes[i].type == 'o':
                NOh -= 1
                Probs = (NDh + NLh)*[self.phi] + NOh*[1]
            Psum = np.cumsum(Probs)
            Psum = Psum/Psum[-1]
            u = np.random.rand(Ninf_cons)
            for j in np.arange(Ninf_cons):
                k = np.sum(u[j]>Psum)
                if k < NDh:
                    infectious_contact = drivers[drivers != i][k]
                else:
                    if k < NDh + NLh:
                        infectious_contact = loaders[loaders != i][k-NDh]
                    else:
                        if k < NDh + NLh + NOh:
                            infectious_contact = office[office != i][k-NDh-NLh]
                    
                if self.nodes[infectious_contact].state == 0:
                    new_infections.append(infectious_contact)
                    self.nodes[infectious_contact].infect_node(day_no,i)
            
            #check for pair infection
            if self.pairs:
                u1 = int(np.random.binomial(1,self.alpha))
                if u1:
                    #find partner
                    for p in d_pairs:
                        if i in p:
                            if p[0] == i:
                                self.nodes[p[1]].infect_node(day_no,i)
                                new_infections.append(p[1])
                            else:
                                self.nodes[p[0]].infect_node(day_no,i)
                                new_infections.append(p[0])
            
        return new_infections
        
    def draw_driver_pairs(self, drivers_available, wp_occ):
        pairs_needed = int(np.round(wp_occ*self.ND/2))   #number of pairs needed for shift
        dpairs = []
        if self.fixed_pairs:     #choose fixed pairings where possible
            pairs_available = []
            drivers_left = drivers_available.copy()   
            for i in np.arange(len(drivers_available)):   #find all fixed pairs of drivers not isolating
                if drivers_available[i]%2 == 0:
                    if i < len(drivers_available)-1 and drivers_available[i+1] == drivers_available[i]+1:
                        pairs_available.append(drivers_available[i]/2)
                        drivers_left.remove(drivers_available[i])
                        drivers_left.remove(drivers_available[i+1])
            
            if len(pairs_available) < pairs_needed:  #if not enough fixed pairs available
                for i in pairs_available:   #all available fixed pairs are used
                    dpairs.append([int(2*i),int(2*i+1)])
                pairs_needed -= len(dpairs)
                if len(drivers_left) < 2*pairs_needed:  #if there are enough drivers, rest are assigned to pairs randomly
                    Nsel = 2*int(len(drivers_left)/2)
                    dsel = np.random.choice(drivers_left, Nsel, replace = False)
                else:   #otherwise as many as possible are assigned randomly
                    dsel = np.random.choice(drivers_left,2*pairs_needed, replace = False)
                for i in np.arange(int(len(dsel)/2)):
                    dpairs.append([int(dsel[2*i]),int(dsel[2*i+1])])
            else:  #if enough fixed pairs available
                if len(pairs_available) == pairs_needed:  #if exactly enough all are assigned
                    Pdraw = pairs_available
                else:    #otherwise enough pairs are chosen at random
                    Pdraw = np.random.choice(pairs_available,pairs_needed, replace = False)
                for i in Pdraw:
                    dpairs.append([int(2*i),int(2*i+1)])   
        else:   #otherwise pairs are chosen at random each day
            if len(drivers_available) < 2*pairs_needed:  
                dsel = np.random.choice(drivers_available, 2*int(len(drivers_available)/2), replace = False)
            else:
                dsel = np.random.choice(drivers_available, 2*pairs_needed, replace = False)
            for i in np.arange(int(len(dsel)/2)):
                dpairs.append([int(dsel[2*i]),int(dsel[2*i+1])])
        
        return dpairs
                                        
    def generate_contacts(self):
        pass
    
    
    def count_states(self):
        SEPIR = np.zeros((5,3))
        in_isolation = np.zeros(3)
        for i in np.arange(self.ND):
            SEPIR[self.nodes[i].state,0] += 1
            if self.nodes[i].isolating:
                in_isolation[0] += 1
        for i in np.arange(self.NL):
            SEPIR[self.nodes[self.ND + i].state,1] += 1
            if self.nodes[self.ND + i].isolating:
                in_isolation[1] += 1
        for i in np.arange(self.NO):
            SEPIR[self.nodes[self.ND + self.NL + i].state,2] += 1
            if self.nodes[self.ND + self.NL + i].isolating:
                in_isolation[2] += 1

        return SEPIR, in_isolation