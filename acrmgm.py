# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 18:28:17 2021
"""
import numpy as np
from scipy.special import loggamma, gamma
from mpmath import hyp3f2
from scipy.stats import norm
from scipy.stats import t as student
import random
from scipy.stats import multivariate_normal

def initialize(n):
    V1 = np.random.choice([0,1], n[0])
    V2 = np.random.choice([0,1], n[1])
    
    c1 = np.zeros( n[0] ).astype(int)
    c2 = np.zeros( n[1] ).astype(int)
    c1[V1==1] = np.random.choice( range( int( np.log( sum(V1==1)))), sum(V1==1))
    c2[V2==1] = np.random.choice( range( int( np.log( sum(V2==1)))), sum(V2==1))
    maxlab = max(max(c1),max(c2)) + 1
    c1[V1==0] = np.random.choice(range(maxlab, maxlab + int( np.log(len(V1) - sum(V1==1))) ),\
                             sum(V1==0)) 
    maxlab = max(c1) + 1
    c2[V2==0] = np.random.choice(range(maxlab, maxlab + int( np.log(len(V2) - sum(V2==1))) ),\
                             sum(V2==0))
    """
    to check that the initialization is in the support 
    verify that common clusters correspond to v=1 using
    V1[np.nonzero(np.in1d(c1, c2))[0]]
    V2[np.nonzero(np.in1d(c2, c1))[0]]
    that have to returns vectors of 1.
    """
    theta1 = np.zeros( n[0] )
    theta2 = np.zeros( n[1] )

    tau1 = np.ones( n[0] )
    tau2 = np.ones( n[1] )
    
    return([c1, c2, V1, V2, theta1, theta2, tau1, tau2])

def samplepseudotie(nc1, nc2, mean1, mean2, var1, var2, \
                     alpha1, alpha2, beta1, beta2, mu1, mu2, lambda1, lambda2,\
                         rho, M):
    a = alpha2 + nc2 / 2
    b = beta2 + 1/2 * nc2 * var2 +\
            nc2 * lambda2 / (2 * (lambda2 + nc2)) * (mean2 - mu2) ** 2
    newtau2 = np.random.gamma(a, 1/b)
        
    m = nc2 / (nc2 + lambda2) * mean2 + \
        lambda2 / (nc2 + lambda2) * mu2
    t = nc2 * newtau2 + lambda2 * newtau2
    newtheta2 = np.random.normal(m, t ** (-1/2))
    
    a = alpha1 + nc1 / 2
    b = beta1 + 1/2 * nc1 * var1 +\
        nc1 * lambda1 / (2 * (lambda1 + nc1 * (1- rho**2))) * (mean1 - mu1) ** 2
    candidates = np.random.gamma(a, 1/b, size = M)
    #
    pcandidates = np.exp(- (mean1 - mu1) * (newtheta2 - mu2) * rho * \
                      nc1 * (lambda1 * lambda2 * newtau2) ** (1/2) /\
                          (lambda1 + nc1 * (1-rho**2)) * (candidates ** (1/2)))
    pcandidates = pcandidates / np.sum(pcandidates)
    newtau1 = np.random.choice(candidates, p = pcandidates)
    
    if (rho!=1):
        m = nc1 /  (nc1  + lambda1 / (1 - rho**2) ) * mean1 + \
        (lambda1 / (1 - rho**2)) / (nc1 + lambda1 / (1 - rho**2) ) *\
        (mu1 + rho * ((lambda2 * newtau2) / (lambda1 * newtau1))**(1/2) *\
        (newtheta2 - mu2))
        t = nc1 * newtau1 + lambda1 * newtau1 / (1 - rho**2) 
    else:
        m = nc1 /  (nc1  + lambda1 ) * mean1 + \
        lambda1 / (nc1 + lambda1 ) *\
        (mu1 + rho * ((lambda2 * newtau2) / (lambda1 * newtau1))**(1/2) *\
        (newtheta2 - mu2))
        t = nc1 * newtau1 + lambda1 * newtau1
    newtheta1 = np.random.normal(m, t ** (-1/2))
        
    return(newtheta1, newtau1, newtheta2, newtau2)

def sampletie(nc, mean, var, alpha, beta, mu, lamb):
    a = alpha + nc / 2
    b = beta + 1/2 * nc * var +\
            nc * lamb / (2 * (lamb + nc)) * (mean - mu) ** 2
    newtau = np.random.gamma(a, 1/b)
        
    m = nc / (nc + lamb) * mean + \
        lamb / (nc + lamb) * mu
    t = nc * newtau + lamb * newtau
    newtheta = np.random.normal(m, t ** (-1/2))
    
    return(newtheta, newtau)

def samplepsuedotie_conditional(nc1, mean1, var1, theta2, tau2,\
                     alpha1, beta1, mu1, mu2, lambda1, lambda2,\
                     rho, M):
    a = alpha1 + nc1 / 2
    b = beta1 + 1/2 * nc1 * var1 +\
        nc1 * lambda1 / (2 * (lambda1 + nc1 * (1- rho**2))) * (mean1 - mu1) ** 2
    candidates = np.random.gamma(a, 1/b, size = M)
    pcandidates = np.exp(-(mean1 - mu1) * (theta2 - mu2) * rho * \
                      nc1 * (lambda1 * lambda2 * tau2) ** (1/2) /\
                          (lambda1 + nc1 * (1-rho**2)) * candidates ** (1/2))
    pcandidates = pcandidates / np.sum(pcandidates)
    newtau1 = np.random.choice(candidates, p = pcandidates)
    
    #m = nc1 /  (nc1  + lambda1 / (1 - rho**2) ) * mean1 + \
    #    (lambda1 / (1 - rho**2)) / (nc1 + lambda1 / (1 - rho**2) ) *\
    #    (mu1 + rho * ((lambda2 * tau2) / (lambda1 * newtau1))**(1/2) *\
    #    (theta2 - mu2))
    #t = nc1 * newtau1 + lambda1 * newtau1 / (1 - rho**2) 
    if (rho!=1):
        m = nc1 /  (nc1  + lambda1 / (1 - rho**2) ) * mean1 + \
        (lambda1 / (1 - rho**2)) / (nc1 + lambda1 / (1 - rho**2) ) *\
        (mu1 + rho * ((lambda2 * tau2) / (lambda1 * newtau1))**(1/2) *\
        (theta2 - mu2))
        t = nc1 * newtau1 + lambda1 * newtau1 / (1 - rho**2) 
    else:
        m = nc1 /  (nc1  + lambda1 ) * mean1 + \
        lambda1 / (nc1 + lambda1 ) *\
        (mu1 + rho * ((lambda2 * tau2) / (lambda1 * newtau1))**(1/2) *\
        (theta2 - mu2))
        t = nc1 * newtau1 + lambda1 * newtau1
    newtheta1 = np.random.normal(m, t ** (-1/2))
    
    return(newtheta1, newtau1)
#p. 17 Bernoulli
def sampleVstar(n, nbar, nc, mbar, theta, z):
    p = np.zeros(2)
    p[0] = np.exp ( loggamma(theta + n[0] - nbar - nc + n[1]) - \
                    loggamma(theta + n[0] - nbar - nc) ) * \
           (z) * \
           hyp3f2 (theta + n[1] - mbar - theta  * z + n[0] - nbar - nc,\
                   n[0], \
                   n[1], \
                   theta + n[1] - mbar + n[0],\
                   theta + n[0] + n[1] - nbar - nc , 1)
    p[1] = np.exp ( loggamma(theta + n[0] - nbar + n[1]) - \
                    loggamma(theta + n[0] - nbar ) ) * \
            (1-z) *\
           hyp3f2 (theta + n[1] - mbar - theta  * z + n[0] - nbar,\
                   n[0], \
                   n[1], \
                   theta + n[1] - mbar + n[0],\
                   theta + n[0] + n[1] - nbar, 1)
               
    p = p / np.sum(p)
    return(np.random.choice([0,1],p=p))

def samplelabel1(i, V1, c1, V2, c2, theta1, theta2, tau1, tau2, \
                 alpha1, alpha2, beta1, beta2, mu1, mu2, lambda1, lambda2,
                 rho, z, theta, simtau1p, X1):
    
    dropc = np.delete(c1, i)
    
    if V1[i] == 1:
        p = np.zeros( (len( np.union1d(dropc [np.delete(V1, i) == 1],\
                                           c2[V2 == 1]) ) + 1, 2) )
        p[:, 0] = np.append( np.unique( dropc [np.delete(V1, i) == 1] ), \
                             np.append( np.setdiff1d(c2[V2 == 1],\
                                  dropc [np.delete(V1, i) == 1]), -1) )
            
        for j in range( len(p) - 1 ):
            nc = np.count_nonzero( dropc == p[j, 0] ) +\
                 np.count_nonzero( c2 == p[j, 0] )
            if np.count_nonzero( dropc == p[j, 0] ) > 0:
                thetac = theta1[c1 == p[j, 0]][0]
                tauc = tau1[c1 == p[j, 0]][0]
                p[j, 1] = nc * norm.pdf(X1[i], loc = thetac, scale = tauc**(-1/2))
            else:
                thetac = theta2[c2 == p[j, 0]][0]
                tauc = tau2[c2 == p[j, 0]][0]
                mutemp = mu1 + (thetac - mu2) * rho *\
                    (lambda2 * tauc / lambda1) ** (1/2) *\
                    simtau1p ** (-1/2)
                tautemp = simtau1p * lambda1 * (1 + lambda1 - rho**2) ** (-1)
                p[j, 1] = nc * np.mean(norm.pdf(X1[i], loc = mutemp, \
                                                scale = tautemp**(-1/2)))
        #to optimaze the code: this are always the same at each iter except z
        p[len(p) - 1, 1] = (1 - z) * theta * student.pdf(X1[i], df = 2 * alpha1,\
                                         loc = mu1, \
            scale = (beta1 * (lambda1 + 1) / (alpha1 * lambda1))**(1/2))
    else:
        p = np.zeros( ( len( np.unique( dropc[np.delete(V1, i) == 0] ) ) + 1, 2 ) )
        p[:, 0] = np.append( np.unique( dropc[np.delete(V1, i) == 0] ), -1)
        
        for j in range( len(p) - 1 ):
            nc = np.count_nonzero( dropc == p[j, 0] ) 
            thetac = theta1[c1 == p[j, 0]][0]
            tauc = tau1[c1 == p[j, 0]][0]
            p[j, 1] = nc * norm.pdf(X1[i], loc = thetac, scale = tauc**(-1/2))
        #to optimaze the code: this are always the same at each iter except z
        p[len(p) - 1, 1] = z * theta * student.pdf(X1[i], df = 2 * alpha1,\
                                         loc = mu1, \
            scale = (beta1 * (lambda1 + 1) / (alpha1 * lambda1))**(1/2)) 
        
    p[:,1] = p[:,1] / np.sum(p[:,1])
    newc = np.random.choice(p[:,0], p = p[:,1])
    return(newc)

def samplelabel2(i, V1, c1, V2, c2, theta1, theta2, tau1, tau2, \
                 alpha1, alpha2, beta1, beta2, mu1, mu2, lambda1, lambda2,
                 rho, z, theta, simtau2p, X2):        
    dropc = np.delete(c2, i)
        
    if V2[i] == 1:
        p = np.zeros( (len( np.union1d(dropc [np.delete(V2, i) == 1],\
                                       c1[V1 == 1]) ) + 1, 2) )
        p[:, 0] = np.append( np.unique( dropc [np.delete(V2, i) == 1] ), \
                             np.append( np.setdiff1d(c1[V1 == 1],\
                                  dropc [np.delete(V2, i) == 1]), -1) )
            
        for j in range( len(p) - 1 ):
            nc = np.count_nonzero( dropc == p[j, 0] ) +\
                 np.count_nonzero( c1 == p[j, 0] )
            if np.count_nonzero( dropc == p[j, 0] ) > 0:
                thetac = theta2[c2 == p[j, 0]][0]
                tauc = tau2[c2 == p[j, 0]][0]
                p[j, 1] = nc * norm.pdf(X2[i], loc = thetac, scale = tauc**(-1/2))
            else:
                thetac = theta1[c1 == p[j, 0]][0]
                tauc = tau1[c1 == p[j, 0]][0]
                mutemp = mu2 + (thetac - mu1) * rho *\
                    (lambda1 * tauc / lambda2) ** (1/2) *\
                    simtau2p ** (-1/2)
                tautemp = simtau2p * lambda2 * (1 + lambda2 - rho**2) ** (-1)
                p[j, 1] = nc * np.mean(norm.pdf(X2[i], loc = mutemp, scale = tautemp**(-1/2)))
        p[len(p) - 1, 1] = (1-z) * theta * student.pdf(X2[i], df = 2 * alpha2,\
                                         loc = mu2, \
            scale = (beta2 * (lambda2 + 1) / (alpha2 * lambda2))**(1/2))
    else:
        p = np.zeros( ( len( np.unique( dropc [np.delete(V2, i) == 0] ) ) + 1, 2 ) )
        p[:, 0] = np.append( np.unique( dropc [np.delete(V2, i) == 0] ), -1)
        
        for j in range( len(p) - 1 ):
            nc = np.count_nonzero( dropc == p[j, 0] ) 
            thetac = theta2[c2 == p[j, 0]][0]
            tauc = tau2[c2 == p[j, 0]][0]
            p[j, 1] = nc * norm.pdf(X2[i], loc = thetac, scale = tauc**(-1/2))
                
        p[len(p) - 1, 1] = z * theta * student.pdf(X2[i], df = 2 * alpha2,\
                                         loc = mu2, \
            scale = (beta2 * (lambda2 + 1) / (alpha2 * lambda2))**(1/2)) 
    
    p[:,1] = p[:,1] / np.sum(p[:,1])
    
    newc = np.random.choice(p[:,0], p = p[:,1])
    return(newc)

def samplez(z, n, c1, c2, V1, V2, theta, alphaz, betaz, eps):
    nbar = np.count_nonzero( V1==0 )
    mbar = np.count_nonzero( V2==0 )
    H1 = len( np.unique( c1[V1==0] ) )
    H2 = len( np.unique( c2[V2==0] ) )
    H = len( np.setdiff1d( c1[V1==1], c2[V2==1] ) ) + \
        len( np.setdiff1d( c2[V2==1], c1[V1==1] ) )
    
    newz = np.random.uniform(z-eps, z+eps)
    pnew = hyp3f2 (theta + n[1] - mbar - theta  * newz + n[0] - nbar,\
                   n[0], \
                   n[1], \
                   theta + n[1] - mbar + n[0],\
                   theta + n[0] + n[1] - nbar, 1) *\
                    newz ** (H1 + H2 + alphaz) * (1 - newz) ** (H + betaz) 
    pold = hyp3f2 (theta + n[1] - mbar - theta  * z + n[0] - nbar,\
                   n[0], \
                   n[1], \
                   theta + n[1] - mbar + n[0],\
                   theta + n[0] + n[1] - nbar, 1) *\
                    z ** (H1 + H2 + alphaz) * (1 - z) ** (H + betaz)
    if(newz < 1 and newz>0 and random.random() < pnew / pold):
        z = newz
    return(z)

def sampletheta(theta, n, c1, c2, V1, V2, z, alphatheta, betatheta):
    nbar = np.count_nonzero( V1==0 )
    mbar = np.count_nonzero( V2==0 )
    K = len( np.union1d( c1, c2 ) )
    
    newtheta = np.random.gamma(alphatheta, 1/ betatheta)
    pnew = np.exp ( loggamma(newtheta + n[0] - nbar) - \
                    loggamma(newtheta + n[0] + n[1] - nbar) ) * \
           hyp3f2 (newtheta + n[1] - mbar - newtheta  * z + n[0] - nbar,\
                   n[0], \
                   n[1], \
                   newtheta + n[1] - mbar + n[0], \
                   newtheta + n[0] + n[1] - nbar, 1) *\
               newtheta **  K
               
    pold = np.exp ( loggamma(theta + n[0] - nbar) - \
                    loggamma(theta + n[0] + n[1] - nbar) ) * \
           hyp3f2 (theta + n[1] - mbar - theta  * z + n[0] - nbar,\
                   n[0], \
                   n[1], \
                   theta + n[1] - mbar + n[0], \
                   theta + n[0] + n[1] - nbar, 1) *\
               theta **  K
    if(pold == 0 and pnew > 0):
        theta = newtheta
    elif(random.random() < pnew / pold):
        theta = newtheta
    return(theta)
        
def samplerho(rho, c1, c2, theta1, theta2, tau1, tau2, lambda1, lambda2, mu1, mu2, eps):
    newrho = random.uniform(rho-eps, rho+eps)
    if (newrho<1 and newrho>-1):
    #if (newrho<-0.5 and newrho>-1):
        pold = 1 
        pnew = 1    
        #note that if len(np.intersect1d(c1, c2))==0,
        #rho is correctly sampled from the prior
        for j in np.intersect1d(c1, c2):
            thetastar1 = theta1[c1 == j][0]
            taustar1 = tau1[c1 == j][0]
            thetastar2 = theta2[c2 == j][0]
            taustar2 = tau2[c2 == j][0]
            
            cov = (lambda1 * taustar1 * lambda2 * taustar2)**(-1/2) * rho
            newcov = (lambda1 * taustar1 * lambda2 * taustar2)**(-1/2) * newrho
            pnew = pnew * \
                multivariate_normal.pdf([thetastar1, thetastar2],\
                                        mean = [mu1, mu2], \
                                    cov = [[(lambda1 * taustar1)**(-1), newcov], \
                                          [newcov, (lambda2 * taustar2)**(-1)]])
            pold = pold * \
                multivariate_normal.pdf([thetastar1, thetastar2],\
                                        mean = [mu1, mu2], \
                                    cov = [[(lambda1 * taustar1)**(-1), cov], \
                                          [cov, (lambda2 * taustar2)**(-1)]])
        if(random.random() < pnew / pold):
            rho = newrho
    return(rho)

def predweights(V1, V2, theta, z, M):
    nbar = np.count_nonzero( V1==0 )
    mbar = np.count_nonzero( V2==0 )
    
    w1 = np.random.beta(nbar + theta*z, len(V1) - nbar + theta, size = M)
    w2 = np.random.beta(mbar + theta*z, len(V2) - mbar + theta, size = M)
    
    v1newtemp = [np.mean(gamma(theta)**2 / gamma(theta*(1-z)) / gamma(theta*(1-z))\
                 * w1 / (1-w1*w2)**(theta*(1-z))),
                 np.mean(gamma(theta)**2 / gamma(theta*(1-z)) / gamma(theta*(1-z))\
                 * (1-w1) / (1-w1*w2)**(theta*(1-z)))]
    v2newtemp = [np.mean(gamma(theta)**2 / gamma(theta*(1-z)) / gamma(theta*(1-z))\
                 * w2 / (1-w1*w2)**(theta*(1-z))),
                 np.mean(gamma(theta)**2 / gamma(theta*(1-z)) / gamma(theta*(1-z))\
                 * (1-w2) / (1-w1*w2)**(theta*(1-z)))]
    v1new = v1newtemp / np.sum(v1newtemp)
    v2new = v2newtemp / np.sum(v2newtemp)
    
    return(v1new, v2new)

def predictive(grid, c1, c2, V1, V2, theta1, theta2, tau1, tau2, \
               alpha1, alpha2, beta1, beta2, lambda1, lambda2, mu1, mu2,
               theta, z, rho, simtau1p):
    dens0 = np.ones(len(grid))
    dens1 = np.ones(len(grid))
    for i in range(len(grid)):
        #v=0
        np0 = len(theta1[V1==0])
        dens0[i] = (np.sum(norm.pdf(grid[i], loc = theta1[V1==0], \
                                      scale = tau1[V1==0]**(-1/2))) +
             theta * z * student.pdf(grid[i], df = 2 * alpha1,\
                                             loc = mu1, \
                scale = (beta1 * (lambda1 + 1) / (alpha1 * lambda1))**(1/2)))\
                 / (np0 + theta * z) 
        #v=1
        p = np.zeros( (len( np.union1d( c1[V1 == 1],\
                                        c2[V2 == 1]) ) + 1, 2) )
        p[:, 0] = np.append( np.unique( c1[V1 == 1] ), \
                             np.append( np.setdiff1d(c2[V2 == 1],\
                                      c1[V1 == 1]), -1) )
                
        for j in range( len(p) - 1 ):
            nc = np.count_nonzero( c1 == p[j, 0] ) +\
                 np.count_nonzero( c2 == p[j, 0] )
            if np.count_nonzero( c1 == p[j, 0] ) > 0:
                thetac = theta1[c1 == p[j, 0]][0]
                tauc = tau1[c1 == p[j, 0]][0]
                p[j, 1] = nc * norm.pdf(grid[i], loc = thetac, scale = tauc**(-1/2))
            else:
                thetac = theta2[c2 == p[j, 0]][0]
                tauc = tau2[c2 == p[j, 0]][0]
                mutemp = mu1 + (thetac - mu2) * rho *\
                    (lambda2 * tauc / lambda1) ** (1/2) *\
                    simtau1p ** (-1/2)
                tautemp = simtau1p * lambda1 * (1 + lambda1 - rho**2) ** (-1)
                p[j, 1] = nc * np.mean(norm.pdf(grid[i], loc = mutemp, \
                                                scale = tautemp**(-1/2)))
        #to optimaze the code: this are always the same at each iter except z
        p[len(p) - 1, 1] = (1 - z) * theta * student.pdf(grid[i], df = 2 * alpha1,\
                                         loc = mu1, \
            scale = (beta1 * (lambda1 + 1) / (alpha1 * lambda1))**(1/2))
                
        np1 = len(theta1[V1==1]) + len(theta2[V2==1])
        
        dens1[i] = np.sum(p[:,1])/(np1 + theta * (1-z)) 
    return(dens0, dens1)

f1_1 = lambda w1, w2, nbar, mbar, n0, n1, theta, z: \
    w1 ** (nbar + theta * z - 1) * w2 ** (mbar + theta * z - 1) * \
    (1 - w1) ** (n0 - nbar + theta ) * (1 - w2) ** (n1 - mbar + theta - 1)/\
    (1 - w1 * w2) ** (theta + theta * z)    

f1_0 = lambda w1, w2, nbar, mbar, n0, n1, theta, z: \
    w1 ** (nbar + theta * z) * w2 ** (mbar + theta * z - 1) * \
    (1 - w1) ** (n0 - nbar + theta - 1) * (1 - w2) ** (n1 - mbar + theta - 1)/\
    (1 - w1 * w2) ** (theta + theta * z)  
    
f2_1 = lambda w1, w2, nbar, mbar, n0, n1, theta, z: \
    w1 ** (nbar + theta * z - 1) * w2 ** (mbar + theta * z - 1) * \
    (1 - w1) ** (n0 - nbar + theta - 1) * (1 - w2) ** (n1 - mbar + theta)/\
    (1 - w1 * w2) ** (theta + theta * z)    

f2_0 = lambda w1, w2, nbar, mbar, n0, n1, theta, z: \
    w1 ** (nbar + theta * z - 1) * w2 ** (mbar + theta * z) * \
    (1 - w1) ** (n0 - nbar + theta - 1) * (1 - w2) ** (n1 - mbar + theta - 1)/\
    (1 - w1 * w2) ** (theta + theta * z)