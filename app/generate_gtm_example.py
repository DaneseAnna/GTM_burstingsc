if __name__ == "__main__":
    import pandas as pd
    import numpy as np
    import random
    from scipy.stats import loguniform, lognorm
    from scipy import stats
    import time
    from os import listdir
    from os.path import isfile, join

    # from multiprocessing import Pool
    from joblib import Parallel, delayed

    import sys
    import statisUtils as utils
    import ABCSampler as abcSampler

    input_file = str(sys.argv[1])
    output_path = str(sys.argv[2])



    def lognpdf(x, mu, sigma):
        shape  = sigma
        loc    = 0
        scale  = np.exp(mu)
        return lognorm.pdf(x, shape, loc, scale)


    df = pd.read_csv(input_file, header=None, index_col=[0]) # gene x cells
    existing_genes = [(f.split('.')[0]).split('_')[1] for f in listdir(output_path) if isfile(join(output_path, f))]


    def parallel_process(args):
        gene_name, row = args

        if gene_name in existing_genes:
            print(f"Skip gene: {gene_name}")
            return 

        start_processing = time.time()
        print(f"Start Processing gene: {gene_name}")

        data = np.array(row)
        
        gene = data[data>=0]
        k = int(0.95*len(gene))
        gene = gene[np.argpartition(gene, k)[:k]]
        gene = np.sort(gene)
        
        statis_data = utils.statisData(gene)
        
        rho = lambda s: np.sqrt(
                    np.sum(
                        np.log(statis_data/(s))**2
                    )
        )
        
        f = lambda k: utils.statisGTM(k,4)
        N = 1000 # the number of particles
        T = 5 # number of rounds

        epsilon = 1# a sequence of discrepancy acceptance thresholds

        prior = lambda: np.array([
                    5 * np.random.uniform(), # kon ~ U[0,5]
                    10**random.uniform(-1,1), # ron ~ logU[-1,1] base 10
                    5 * np.random.uniform(), # koff ~ U[0,5]
                    10**random.uniform(-1,1), # roff ~ logU[-1,1] base 10 
                    100 * np.random.uniform(), # rsyn ~ U[0,100] 
                    1]) # rdeg=1

        proposal_sigma = 0.2
        
        proposal = lambda x: np.random.lognormal(x,proposal_sigma)
        
        proposal_pdf = lambda kon_post,kon_prior,ron_post,ron_prior,koff_post,koff_prior,roff_post,roff_prior,mu_post,mu_prior: \
        lognpdf(mu_post,np.log(mu_prior),proposal_sigma) * lognpdf(kon_post,np.log(kon_prior),proposal_sigma) \
        * lognpdf(ron_post,np.log(ron_prior),proposal_sigma) * lognpdf(koff_post,np.log(koff_prior),proposal_sigma) \
        * lognpdf(roff_post,np.log(roff_prior),proposal_sigma)

        result, flag = abcSampler.ABCSMCSampler(N,prior,f,rho,epsilon,T,proposal,proposal_pdf,gene_name)

        np.save(f'{output_path}/gene_{gene_name}.npy', np.array([gene_name, gene, result], dtype=object))

        end_processing = time.time()
        total_processing = end_processing - start_processing
        print(f"Finish Processing gene: {gene_name}, time= {total_processing}")


    # run script in parallel
    res = Parallel(
        n_jobs=5
    )(
        delayed(parallel_process)((index, row)) for index, row in df.iterrows()
    )

