
'''
work:
1. run higher terms in DGX1
2. learn second order monte carlo
3. use "if" to realize the boundary condition in the paper

'''
import numpy as np
import matplotlib.pyplot as plt



#all parameters large than 0
min_v = 0.2
max_v = 2
min_price = 1 #since the uniform mesh, then dv = dx
max_price = 40
rho = 0.5
theta = 1  #mean min vol and max vol
kappa = 2  # the degree of mean reversion
r = 0.2#0.02 #risk-free rate
sigma = 1 #volvol

K = 20    #strike price 


T = 1
L = 1
Ns = 1000
Nv = 500 #space or stock price interval(less than sqrt(Nt))
Nt = 400000 #time

ds = (max_price - min_price)/Ns #1
dv = (max_v - min_v)/Nv #1
dt = T/Nt



u = np.zeros((Nt+1, Ns+1, Nv+1))

if dt*(r*(dv**2)*(ds**2) + max_v*(sigma**2)*(ds**2) + max_v*(max_price**2)*(dv**2)) <= (dv**2) * (ds**2):
    #dt <= ((dv**2) * (ds**2))/(r*(dv**2)*(ds**2) + max_v*(sigma**2)*(ds**2) + max_v*(max_price**2)*(dv**2)):
    #initial condition applied
    for i in range(0, Ns+1):
        for j in range(0, Nv+1):
            u[0,i,j] = max((min_price+(i*ds)-K),0)  
            

            
    for t in range(0, Nt):
        
        #old boundary condition
        #top
        u[t+1,0,:] = 0 
        #left
        u[t+1,0,:] = 0
        #botton
        u[t+1,:,Nv] = 0
        
        u[t+1,Ns,:] = 0

        
        
        for i in range(1, Ns):
            si = min_price + i*ds

            for j in range(1, Nv):  
                vj = min_v + j*dv
                
            
                '''
                #new boundary condition
                u[t,0,j] = 0 
                u[t,Ns,j] = ds + u[t,Ns-1,j]
                u[t+1,i,0] = u[t,i,0] - (r*si*dt/ds)*(u[t,i+1,0] - u[t,i,0]) - (kappa*theta*dt/dv)*(u[t,i,1] - u[t,i,0]) + dt*r*u[t,i,0]
                u[t,i,Nv] = si
                '''
                
                #iteration
                a1 = 1 - (r*dt) - (dt*vj*(si**2))/(ds**2) - (dt*vj*(sigma**2) / (dv**2))
                a2 = (dt/2)*si*(vj*si/(ds**2) + (r/ds))
                a3 = (dt/2)*si*(vj*si/(ds**2) - (r/ds))
                a4 = (dt/2)*((sigma**2)*vj/(dv**2) + kappa*(theta - vj)/dv)
                a5 = (dt/2)*((sigma**2)*vj/(dv**2) - kappa*(theta - vj)/dv)
                a6 = (dt*rho*sigma*vj*si)/(4*ds*dv)
                u[t+1,i,j] = a1*u[t,i,j] + a2*u[t,i+1,j] + a3*u[t,i-1,j] + a4*u[t,i,j+1] + a5*u[t,i,j-1] + a6*(u[t,i+1,j+1] - u[t,i+1,j-1] - u[t,i-1,j+1] + u[t,i-1,j-1])
                if i == (Ns-1):
                    u[t+1,Ns,j] = ds + u[t+1,Ns-1,j]
                
    row = np.linspace(min_price, max_price, Ns+1)
    columns = np.linspace(min_v, max_v, Nv+1) 
    columns, row= np.meshgrid(columns, row)
    ax = plt.axes(projection='3d')
    ax.contour3D(columns, row, u[-1], 50)
    ax.set_xlabel('violatility (v)')
    ax.set_ylabel('stock price (s)')
    ax.set_zlabel('option price')
    ax.set_title('Heston option price model')

else:
    print('stability condition failed')
