import numpy as np
from scipy import stats
import time


def ABCSMCSampler(N,prior,f,rho,epsilon,T,proposal,proposal_pdf,gene):
    """Sequential Monte Carlo Sampler for approximate Bayesian computaion

    Inputs:
        N - the number of particles
        prior - prior distribution sampler
        f - function that generates simulated data give a parameters set
        rho - discrepancy metric, treated as a function of simulated data only
        epsilon - a sequence of discrepancy acceptance thresholds
        T - number of rounds
        proposal - proposal kernel process, generates new trials
        proposal_pdf - the PDF of the proposal kernel
        gene - gene-name

    Outputs:
        result - the sequence of particles of size N x T
    """

    # initialize
    result, flag = ABCRejectionSampler(N,prior,f,rho,epsilon,T,gene)
    param = {}
    W = (1 / N) * np.ones((N,T)) # sequential sampling
    
    if flag == True: 
        for t in range(2, T + 1):
            temp_dist_array = np.array([d[0]['dist'] for d in result[:,t-2]]) 
            epsilon = np.median(temp_dist_array) # threshold

            # generate new generation of particles
            for i in range(1, N+1):
                r = np.inf
                total_time = 0

                while total_time < 10 and r > epsilon: # rejections steps
                    start = time.time()

                    j = np.random.choice(N, size = 1, replace = True, p = W[:, t-2])[0]
                    param_temp = result[j,t-2]

                    # generate new particle based on proposal kernel
                    param_proposal = proposal(
                        np.log(
                            [
                                param_temp[0]['kon'],
                                param_temp[0]['ron'],
                                param_temp[0]['koff'],
                                param_temp[0]['roff'],
                                param_temp[0]['mu']
                            ]
                        )
                    )
                    
                    param['kon'] = param_proposal[0]
                    param['ron'] = param_proposal[1]
                    param['koff'] = param_proposal[2]
                    param['roff'] = param_proposal[3]
                    param['mu'] = param_proposal[4]
                    param['delta'] = 1

                    static_temp = f(param)

                    # condition to discard/skip negative or 0 values
                    if (len(static_temp[static_temp<0]) > 0) or (np.sum(static_temp) == 0.):
                        static_temp = np.abs(static_temp)
                        
                    r = rho(static_temp)
                    

                    result[i-1,t-1][0] = {}
                    result[i-1,t-1][0]['kon'] = param['kon']
                    result[i-1,t-1][0]['ron'] = param['ron']
                    result[i-1,t-1][0]['koff'] = param['koff']
                    result[i-1,t-1][0]['roff'] = param['roff']
                    result[i-1,t-1][0]['mu'] = param['mu']
                    result[i-1,t-1][0]['delta'] = 1
                    result[i-1,t-1][0]['dist'] = r

                    end = time.time()
                    elapsedTime = end - start
                    total_time = total_time + elapsedTime
                
                if total_time > 10:
                    flag = False
                    break

                # recompute particle weight using optimal backward kernel
                back_K = 0
                for j in range (1 , N + 1):
                    back_K = back_K + W[j-1,t-2] * proposal_pdf(
                        result[i-1,t-1][0]['kon'],
                        result[j-1,t-2][0]['kon'], 
                        result[i-1,t-1][0]['ron'],
                        result[j-1,t-2][0]['ron'], 
                        result[i-1,t-1][0]['koff'],
                        result[j-1,t-2][0]['koff'], 
                        result[i-1,t-1][0]['roff'],
                        result[j-1,t-2][0]['roff'], 
                        result[i-1,t-1][0]['mu'],
                        result[j-1,t-2][0]['mu']
                    )
                W[i-1,t-1] = ((1/5)*(1/5)*(1/30)*1/(4*result[i-1,t-1][0]['ron'] * result[i-1,t-1][0]['roff']))/back_K
            
            if flag == False:
                break

            # resample
            if t < T:
                result_rs  = result[:,t-1]
                W[:,t-1] = W[:,t-1]/sum(W[:,t-1])
                J =  np.random.choice(N, size = N, replace = True, p = W[:, t-1])
                result[:,t-1] = result_rs[J]
            
            # re-set weights
            W[:,t-1] = 1./N

    return (result, flag)   


def ABCRejectionSampler(N,prior,f,rho,epsilon,T,gene):
    """Rejection Sampler for approximate Bayesian computaion

    Inputs:
        N - the number of ABC posterior samples
        prior - function that generates iid samples from the parameter joint
            prior distribution
        f - function that computes statictis given a parameters set
        rho - discrepancy metric, treated as a function of simulated data only
        epsilon - the discrepancy acceptance threshold
        T - number of rounds
        gene - gene-name

    Outputs:
        result - a matrix of ABC prior samples
    """
    
    result = np.empty([N,T,1], dtype=object)
    result0 = []
    total_time = 0.
    flag = True

    size_val = 0 #indicates the accepted index

    while total_time < 5*N and size_val < 5*N: #rejection criteria
        start = time.time()
        param = {}

        # generate trial from the prior
        theta_trial = prior()
        param['kon'] = theta_trial[0]
        param['ron'] = theta_trial[1]
        param['koff'] = theta_trial[2]
        param['roff'] = theta_trial[3]
        param['mu'] = theta_trial[4]
        param['delta'] = theta_trial[5]

        # compute theorical statictis of parameters
        static_theo = f(param)

        # condition to discard/skip negative or 0 values
        if (len(static_theo[static_theo<0]) > 0) or (np.sum(static_theo) == 0.):
            continue      

        dist = rho(static_theo)
        param['dist'] = dist

        # accept or reject
        if dist <= epsilon:
            result_obj = {}
            result_obj['kon'] = param['kon']
            result_obj['ron'] = param['ron']
            result_obj['koff'] = param['koff']
            result_obj['roff'] = param['roff']
            result_obj['mu'] = param['mu']
            result_obj['delta'] = param['delta']
            result_obj['dist'] = param['dist']
            
            result0.append(result_obj)
            size_val += 1  
        
        end = time.time()
        elapsedTime = end - start
        total_time = total_time + elapsedTime

    if total_time > 5*N:
        flag = False
        print(f'Gene {gene}:wrong!!\n')
    else:
        index0 = np.argsort([d['dist'] for d in result0])
        result0 = np.take(result0, index0[:N])
        for index2 in range (1, N + 1):
            result[index2-1,0] = result0[index2-1]

    return (result, flag)
