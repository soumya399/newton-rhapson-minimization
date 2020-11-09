
# program to calculate the minimum pressure required to break the ice by newton-rhapson method
import matplotlib.pyplot as plt
import math
import numpy as np


#function
def model(p,h):
	return (pow(p,3)*(1-pow(beta,2)))+(((0.4*h*pow(beta,2))-(sigma*pow(h,2)/pow(r,2)))*pow(p,2))+(((pow(sigma,2)*pow(h,4))/(3*pow(r,4)))*p)-(pow((sigma*pow(h,2)/(3*pow(r,2))),3))

#derivative of function
def model_prime(p,h):
	return(3*(pow(p,2)*(1-pow(beta,2))))+(2*p*((0.4*h*pow(beta,2))-(sigma*pow(h,2)/pow(r,2))))+((pow(sigma,2)*pow(h,4))/(3*pow(r,4)))

# input parametters for the equation

beta=0.5

r=40

# Tensile strength of the ice (psi)
sigma = 150

# thickness values for the ice
height = [0.6,1.2,1.8,2.4,3.0,3.6,4.2]

# initial guessed value of the pressure
p_guess = 0.1

# relaxation factor
alpha = 1

# defining the tolerance eroor which is acceptable
tolerance = 1e-4


#Setting up various relaxation factors
relax_factors = [0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5]

h_value = 0.6
iters = []


for rf in relax_factors: 
	it = 1 
	p_guess = 0.1
	while(abs(model(p_guess,h_value))>tolerance):
		p_guess = p_guess-((model(p_guess,h_value)/model_prime(p_guess,h_value))*rf)
		it = it+1


	iters.append(it)
	

print(iters)
plt.plot(relax_factors,iters)
plt.xlabel("relaxation factor")
plt.ylabel("number of iterations")
plt.show()



pressure=[]

iter=1
for h in height:
	while(abs(model(p_guess,h))>tolerance):
		p_guess = p_guess-((model(p_guess,h)/model_prime(p_guess,h))*alpha)
		iter = iter+1
	pressure.append(p_guess)

print ("thickness  pressure")
for i in range(0,7):
		
	print(height[i], "   " ,pressure[i])
print("\nnumber of iterations: ")
print(iter)

plt.plot(height,pressure)
plt.xlabel('thickness (in feet)')
plt.ylabel('pressure(in psi')
plt.title('pressure vs thickness')
plt.xlim([0,6])
plt.ylim([0,0.5])
plt.show()



