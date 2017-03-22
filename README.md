## This project is to study the effect of temperature and the strength of random field on the order of a one-dimensional chain of N harmonically coupled particles.
## This Python code was run Python 2.7 using Xcode Version 6.1.1 in mac os x version 10.9.5.
 
Theoretically, higher temperature and strong random field will destroy the order of the system more.

## Algorithm
The initial configuration for the N-particle chain is chosen at random.

Start doing Monte Carlo steps
1. Pick particle-i randomly from N particles.

2. Pick displacement d for particle-i randomly from an array with displacement values from −γ to γ. Here γ depends on temperature, because for higher temperature, the particles have chance to move further. In my code, I use a linear relationship between γ and temperature.

3. Calculate the energy change for particle-i ΔE = (1/2)\*k[x(i)+d−x(i−1)]^2 + (1/2)\*k[x(i+1)−x(i)−d]^2 +RF (i), where RF(i) is the random field on each particle, and we can change the strength of the random field. Here I use periodic boundary condition, that is particle-N is connected to particle-1 by a spring too. (In the code, particle-N corresponds to x(N − 1) and particle-1 corresponds to x(0).)

4. Compare \exp(−ΔE/T) with a random number R (0\<R\<1). If \exp(−ΔE/T) \> R, accept the move and update x(i) = x(i) + d. Otherwise go to step 1.

5. repeat step 1 to 4 for N times.
