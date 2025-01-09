import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mesa
from mesa import Agent, Model
from mesa.time import RandomActivation
import sys
import random
import socket 

#Define the parameter and what might happen for a single individiual.
#The random genetic drift is therefore introduced on the individual level.
class organism(Agent) :
    def __init__(self, unique_id, model):                       # INITIAL PARAMETERS
        super().__init__(unique_id, model)
        self.genotype = ''                                      #genotype will be assigned during the run
        self.genotype_frequencies = model.genotype_frequencies  #genotype_frequencies are stored here. 
        self.mut_prob = model.mut_prob                          # Mutation probability
        self.reproduction_rate = model.reproduction_rate        # Reproduction rate
        self.survival_prob = model.survival_prob                # dictionary of survival probability per genotype
        self.state = 'alive'                                    # State of the organism - dead or alive
        self.num_offspring = 0                                  # number of offsprings per the organism
        self.offspring = None                                   # will assign an agent if the organism reproduces
        self.r_move = model.r_move                              # the radius of possible area for placing an offspring.
                                                                # if this area is full, there won't be reproduction.

    def move_offspring(self):                                   # GETS LOCATION FOR OFFSPRING

        possible_steps = self.model.grid.get_neighborhood(self.pos,moore=True,include_center=False,radius=self.r_move)

        empty_position = (-1, -1)

        #allpos = {}
        emptykeys = []

        for i in range(len(possible_steps)):

            if self.model.grid.is_cell_empty(possible_steps[i]):

                emptykeys.append(i)

                empty_position = possible_steps[i]


        if empty_position != (-1,-1):

            new_key = random.choice(emptykeys)

            new_position = possible_steps[new_key] #allpos[new_key]

            self.offspring_pos = new_position

        else :

            self.offspring_pos = (-1, -1)

    def mutate(self):                        # RANDOM OCCURENCE OF A MUTATION ON WILD-TYPE ALLELE

        R1 = random.random()

        if R1 < self.mut_prob:

            if self.genotype == '++':

                self.genotype = '+m'

            elif self.genotype == '+m':

                self.genotype = 'mm'

    def reproduce(self):                     # REPRODUCTION OF THE ORGANISM
        Rrep = random.random()

        if Rrep < self.reproduction_rate:    # Randomized probability of reproduction   
            R2 = random.random()

            mate_genotype_prob = self.genotype_frequencies

            if R2 < mate_genotype_prob[0]:    # Randomized mate choice -- defined over the mate genotype  

                mate_genotype = '++'

            elif R2 < np.sum(self.genotype_frequencies[0:1]):

                mate_genotype = '+m'

            else:
                mate_genotype = 'mm'

            R3 = random.randint(0,1)
            R4 = random.randint(0,1)

            allele_1 = mate_genotype[R3]       # Randomized segragation of alleles 
            allele_2 = self.genotype[R4]

            off_genotype = allele_1 + allele_2 # Randomized fertilization

            if off_genotype == 'm+':
                off_genotype = '+m'

            self.off_genotype = off_genotype

            self.num_offspring += 1

            self.offspring = organism

            self.move_offspring()

    def survive(self):                       # SURVIVAL OF THE ORGANISM UNDER NATURAL SELECTION

        Rsur = random.random()

        if Rsur > self.survival_prob[self.genotype]: # Randomized survival over the probability of survival

            self.state = 'dead'

            self.model.grid.remove_agent(self) 



    def activity(self):                    # ACTIVITY OF ORGANISM PER GENERATION

        if self.state == 'alive':
            self.mutate()
            self.reproduce()
            self.survive()

    # iterate cell activity #
    def step(self):                        # ITERATION OF AGENT ACTIVITY IN MODEL

        self.activity()

# The model defines the population dynamics
#Val is a dictionary that contains the initial parameters. See below.
# p : The frequency of the dominant allele. q and genotypic frequences of initial population is estimated over p.
# mut_prob : Probability for a dominant WT allele mutating into the recessive mutant allele.
# survival_prob : PRobability of survival (1-selection_coefficient) for each genotype.
# pop_size : Initial population size
# grid : Creates the map of the habitat.
# track_frequencies : Tracks changes in genotypic and allelic frequencies each generation.
#n_alive : Number of alive organisms each generation.
#r_move : Radius of available locations for an offspring. Introduces "resource limit".

class population(Model):

    def __init__(self,Val):
        super().__init__(Val)
        self.genotype_frequencies = [Val['p']*Val['p'], 2*Val['p']*(1-Val['p']), (1-Val['p'])*(1-Val['p'])]
        self.mut_prob = Val['mutation_probability']
        self.reproduction_rate = Val['reproduction_rate']
        self.survival_prob = Val['survival_prob']
        self.pop_size = Val['initial_size']
        self.schedule = RandomActivation(self)
        self.grid = mesa.space.MultiGrid(Val['width'], Val['height'],True)
        self.cntr = 0
        self.running = True
        self.track_frequencies = {'f_pp' : [self.genotype_frequencies[0]],
                                  'f_pq' : [self.genotype_frequencies[1]],
                                  'f_qq' : [self.genotype_frequencies[2]],
                                  'f_p'  : [Val['p']],
                                  'f_q'  : [1-Val['p']]
                                 }
        self.n_alive = self.pop_size
        self.r_move = 1

        num_pp = np.round(self.pop_size*self.genotype_frequencies[0])
        num_pq = np.round(self.pop_size*self.genotype_frequencies[1])
        num_qq = self.pop_size - (num_pp + num_pq)

        self.distribution = [num_pp,num_pq,num_qq]  # will be used to create organisms of different genotypes with respective ratios.

        for i in range(self.pop_size):

            a = organism(i, self)

            self.schedule.add(a)

            x = self.random.randrange(self.grid.width)  

            y =  self.random.randrange(self.grid.height) 

            self.grid.place_agent(a, (x, y))   # places each organism at the start at a random location.

            if i < self.distribution[0]:  # assigns genotypes to each organism
                a.genotype = '++'

            elif i < (self.distribution[0] + self.distribution[1]):
                a.genotype = '+m'

            else:
                a.genotype = 'mm'



        # Data Collection from organism and population
        self.datacollector = mesa.DataCollector(agent_reporters={ "cell id":"unique_id",
                                                             "genotype":"genotype",
                                                             "number of offsprings": "num_offspring"
                                                             },model_reporters={"number of alive cells" : "n_alive",
                                                                               "genotypic frequencies" : "genotype_frequencies"
                                                                               })

    def update_frequencies(self):                 #UPDATE GENOTYPIC FREQUENCIES EACH GENERATION
        n_alive, n_pp, n_pq, n_qq = 0, 0, 0, 0

        for c in self.schedule.agents:

            if c.state == 'alive':

                n_alive +=1

                if c.genotype == '++':

                    n_pp += 1

                elif c.genotype == '+m':

                    n_pq += 1

                elif c.genotype == 'mm':

                    n_qq += 1

        f_pp = n_pp / n_alive
        f_pq = n_pq / n_alive
        f_qq = n_qq / n_alive
        f_q = np.sqrt(f_qq)
        f_p = 1- f_q

        self.track_frequencies['f_pp'].append(f_pp)
        self.track_frequencies['f_pq'].append(f_pq)
        self.track_frequencies['f_qq'].append(f_qq)
        self.track_frequencies['f_q'].append(f_q)
        self.track_frequencies['f_p'].append(f_p)

        self.genotype_frequencies = [f_pp,f_pq,f_qq]

        self.n_alive = n_alive

    def plot_population(self):
        gens = np.arange(0, self.generations+1, 1).tolist()

        title = ('Fitness : f_PP = ' + str(self.survival_prob['++']) + ', f_Pp = ' + str(self.survival_prob['+m']) +
                ', f_pp = ' + str(self.survival_prob['mm']) + '; mutability :' + str(self.mut_prob) +'; reproduction rate :'
                + str(self.reproduction_rate) + '; N_0 =',str(Val['initial_size']))

        plt.figure
        plt.plot(gens,model.track_frequencies['f_p'],'r',label ='p')
        plt.plot(gens,model.track_frequencies['f_q'],'b', label ='q')
        plt.xlabel('Generations')
        plt.ylabel('Allelic Frequencies')
        plt.legend()
        plt.title(title)
        plt.xlim(left = 0)
        plt.ylim([0,1.1])
        plt.grid()
        plt.show()

        plt.figure
        plt.plot(gens,model.track_frequencies['f_pp'],'r',label ='Homozygous dominant')
        plt.plot(gens,model.track_frequencies['f_pq'],'b', label ='Heterozygous')
        plt.plot(gens,model.track_frequencies['f_qq'],'g', label ='Homozygous recessive')
        plt.xlabel('Generations')
        plt.ylabel('Genotypic Frequencies')
        plt.legend()
        plt.title(title)
        plt.xlim(left = 0)
        plt.ylim([0,1.1])
        plt.grid()
        plt.show()


    def step(self):                   # ITERATION OF POPULATION DYNAMICS

        self.update_frequencies()

        self.datacollector.collect(self)

        self.schedule.step()

        self.cntr += 1

        for c in self.schedule.agents :

            if c.offspring :

                if c.offspring_pos != (-1,-1):

                    d = organism(self.pop_size+1,self)

                    d.genotype = c.off_genotype

                    d.pos = c.offspring_pos

                    self.grid.place_agent(d, d.pos)

                    self.schedule.add(d)

                    self.pop_size += 1

                    c.offspring = None

# Finds free port for running visual simulation
def find_free_port():
    with socket.socket() as s:
        s.bind(('', 0))            # Bind to a free port provided by the host.
        return s.getsockname()[1]

# Creates agent portrayal (organism) with respect to genotypes for visual simulation.
def agent_portrayal(agent):
    portrayal = {
        "Shape": "circle",
        "Filled": "true",
    }

    if agent.genotype == '++':
        portrayal["Color"] = "red"
        portrayal["r"] = 1
        portrayal["Layer"] = 1

    if agent.cell_type == '+m':
        portrayal["Color"] = "blue"
        portrayal["r"] = 1
        portrayal["Layer"] = 1

    if agent.cell_type == 'mm':
        portrayal["Color"] = "green"
        portrayal["r"] = 1
        portrayal["Layer"] = 1

    return portrayal

# Visual Simulation
def visualrun(Val):
    grid = mesa.visualization.CanvasGrid(agent_portrayal, Val['width'], Val['height'], 500, 500)
    server = mesa.visualization.ModularServer(population, [grid], "Population Genetics Simulation", {"Val": Val})
    server.port = find_free_port() 
    print(server.port)
    server.launch()

# Run model. key is 'visual' for visual simulation. Anything else for regular simulation.
def run_population(key):
    Val = {}

    Val['p']                    = float(input('Enter the frequency of dominant allele (p) : \n'))

    Val['mutation_probability'] = float(input('Enter the mutability (probably for random mutation into recessive allele) : \n'))

    Val['reproduction_rate']    = float(input('Enter reproduction rate : \n'))

    s_pp                        = float(input('Enter selection coefficient for homozygous dominant genotype :\n'))
    s_pq                        = float(input('Enter selection coefficient for heterozygous genotype :\n'))
    s_qq                        = float(input('Enter selection coefficient for homozygous recessive genotype :\n'))

    Val['survival_prob']        = {'++' : 1-s_pp,
                                   '+m' : 1-s_pq,
                                   'mm' : 1-s_qq
                                  }

    Val['initial_size']         = int(input('Enter initial population size : \n'))

    Val['width']                = int(input('Enter width : '))

    Val['height']               = int(input('Enter height : '))

    n_grid = Val['width']*Val['height']

    while n_grid < Val['initial_size']:

        print('Population size does not fit the map with width ', Val['width'], ' and height ', Val['height'], '. Reenter the values.\n')

        print('Recommendation : Estimate the size of the map to accomodate population growth.')

        Val['width']                = int(input('Enter width : '))

        Val['height']               = int(input('Enter height : '))

        n_grid = Val['width']*Val['height']

    if key == 'visual':
        visualrun(Val)
    else:
        generations = int(input('Enter number of generations : '))

        model = population(Val)

        model.generations = generations

        for i in range(0,generations):
            model.step()
            ga = model.datacollector.get_agent_vars_dataframe()
            gm = model.datacollector.get_model_vars_dataframe()

        gens = np.arange(0, generations+1, 1).tolist()

        plt.figure
        plt.plot(gens,gm['number of alive cells'],'k')
        plt.xlabel('Generations')
        plt.ylabel('Population Size')
        plt.legend()
        plt.title(title)
        plt.xlim(left = 0)
        plt.ylim(bottom = 0)
        plt.grid()
        plt.show()

        return model, ga, gm
