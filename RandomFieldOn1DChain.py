#!/usr/bin/env python

import numpy
from numpy import *
from pylab import *
from scipy import *
from scipy import weave
# weave was the only module never ported to Python 3.x. Remember to change Edit Scheme/Info/Excutable.
# For Python 3.5, in the terminal, open /usr/local/bin
# For Python 2.7, in the terminal, open /System/Library/Frameworks/Python.framework/Versions/2.7/bin/

from scipy.weave import converters



class oneDChain:
    
    """Class that describes a one-dimensional particles connected by springs"""
        
        
        
    def __init__(self, N = 50, temperature = 5.0,k = 1.0, strength = 0.1):
            
            #numpy.random.seed(222)
            self.N = N # number of particles
            self.k = k # spring constant
            self.temperature = temperature
            self.equilibrium = zeros(self.N, float)
            self.strength = strength
                
            self.gamma=0.5 * self.temperature # the maximum possible displacement dependent on temperature
            self.dis=1.0
                
            increment = 0.5 * self.k *self.dis**2/self.temperature
            numbers = int(int(2*self.gamma/self.dis) * (int(2*self.gamma/self.dis)-1)/2)
            self.w = zeros(numbers, float)
            for i in range(numbers):
                    self.w[i] = exp(-increment * i) # store Boltzmann weights
                
            self.aveDisField = zeros(int(self.N/2), float)
        
        
    def displacement(self): # N possible displacement for each particle
            
            return arange(-self.gamma, self.gamma, 2 * self.gamma/self.N) # displacement can be positive or negative
        
        
    def reset(self):
            
            self.monteCarloSteps = 0
            self.acceptedMoves = 0
            self.energyArray = array([], float)
            self.count = zeros(int(self.N/2),int)
            self.aveDisField = zeros(int(self.N/2), float)
            self.state = numpy.random.random(self.N)
        
    def resetT(self):
            self.gamma=0.5 * self.temperature # the maximum possible displacement dependent on temperature
            self.dis=1.0
        
            increment = 0.5 * self.k *self.dis**2/self.temperature
            numbers = int(int(2*self.gamma/self.dis) * (int(2*self.gamma/self.dis)-1)/2)
            self.w = zeros(numbers, float)
            for i in range(numbers):
                    self.w[i] = exp(-increment * i) # store Boltzmann weights
        
        
    def weaveMonteCarloStep(self):
            
            N = self.N
            k = self.k
            displacement = self.displacement()
            dis = self.dis
            w = self.w
                
            state = self.state
            acceptedMoves = array([self.acceptedMoves], int)
                
                
            length = len(self.displacement())
                
            randomPositions = numpy.random.random(N)*self.N # the particles are chosen at random positions
                
            randomMetro = numpy.random.random(N) # random number to do the Metropolis algorithm
                
            randomDis = numpy.random.random(length)*length # randomly pick the displacement for each particle
                
            strength = self.strength
                
            randomPotential = numpy.random.random(N) * strength # strength of the random field
                
                
            code ="""
                    int N1 = N - 1;
                    for (int l = 0; l < N; l++) {
                    
                    int d = displacement(int(randomDis(l)));
                    int i = int(randomPositions(l));
                    float dE = 0.5 * k * pow((state(i) + d -state((i+N1)%N)),2) + 0.5 * k * pow((state(i+1) -state(i) - d),2) +
                    randomPotential(i);
                    if (( w(int(dE/(0.5 * k * pow(dis, 2)))))> randomMetro(l)) {
                    acceptedMoves(0) += 1;
                    state(i) += d;
                    }
                    
                    }
                    
                    """
                
            weave.inline(code, ['N', 'randomDis','displacement','randomPotential','randomPositions','k','state','dis','w', 'randomMetro', 'acceptedMoves'],
                            type_converters=converters.blitz, compiler='gcc')
                    
            energy = 0.5 * self.k *sum(pow(self.state[1:]-self.state[0:-1],2))
                             
            self.state = state
            self.acceptedMoves = acceptedMoves[0]
                             
            self.energyArray = append(self.energyArray, energy)
            self.monteCarloSteps += 1



    def steps(self, number = 100):
        for i in range(number):
            self.weaveMonteCarloStep()
        
        
        # averaged displacement field
    def averageDisplacementField(self):
            
            N = self.N
            displacementField = self.state-self.equilibrium
                
            count = self.count
            aveDisField = self.aveDisField #averaged displacement field
                
            code = """
                    for (int j = 0; j < int(N/2); j++)
                    {
                    for ( int i = 0; i < N; i++)
                    {
                    aveDisField(j) += pow(displacementField(j+i)-displacementField(i),2);
                    count(j)++ ;
                    }
                    for ( int k = 0; k < N; k++)
                    {int kk = j - k;
                    if (kk > 0)
                    {
                    aveDisField(j) += pow(displacementField(k)-displacementField(N-j+k),2);
                    count(j)++;
                    }
                    }
                    }
                    
                    for (int i = 0; i < int(N/2); i++)
                    {
                    aveDisField(i)=aveDisField(i)/count(i);
                    }
                    
                    """
                
            weave.inline(code, ['N', 'aveDisField', 'displacementField','count'],
                             type_converters=converters.blitz, compiler='gcc')
                             
            return aveDisField
        
        
        
    def plot(self):
            x=arange(1,int(self.N/2),1)
                
            figure(1)
            plot(x,self.averageDisplacementField()[1::],label="T = %2f"%(self.temperature))
            xlabel("relative position of the particles ( r )")
            ylabel("displacement field $(u(r)-u(0))^2)")
            legend(loc="best")
                
                
            figure(2)
            plot(self.state,label="T = %2f"%(self.temperature))
            xlabel("position of the particles")
            ylabel("displacement of each particle")
            legend()
                
                
                
            figure(3)
            plot(self.energyArray,label="T = %2f"%(self.temperature))
            xlabel("number of Monte Carlo steps")
            ylabel("energy")
            legend(loc="best")




Chain=oneDChain(N=100, temperature=6,k=0.2, strength = 0.2)

Chain.reset()
Chain.resetT()
Chain.steps(number=10000)
Chain.plot()
print (Chain.acceptedMoves)

Chain=oneDChain(N=100, temperature=20,k=0.2, strength = 0.2)

Chain.reset()
Chain.resetT()
Chain.steps(number=10000)
Chain.plot()
print (Chain.acceptedMoves)

'''
    for i in range(6):
    numpy.random.seed(222)
    Chain.reset()
    Chain.temperature=Chain.temperature + 1.0
    Chain.resetT()
    Chain.steps(number=10000)
    Chain.plot()
    print Chain.acceptedMoves
    '''

show()