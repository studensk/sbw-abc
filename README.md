# sbw-abc

functions.R contains all functions used to run the ABC SMC
  - post_next.sample: runs iterations 2-T of the SMC
       - sample.discrete: sample from discrete distribution based on weights
       - sample.n: takes either a number of draws (for first iteration) or a dataset of parameters (for subsequent iterations), and generates new samples and results.
          - post_theta.sample: generate dataset from a parameter draw
            - post_eps: uses a sample from the prior to filter for trajectory endpoints
          - tkernel.sample: perturb the current set of parameters using a specified kernel
       - tkernel.density: multivariate kernel density function
hysplitrun_2022.R runs the algorithm with the true observed L2 data; outputs code/output/new_results_test.csv
data_generation.R generates a binary L2 (presence/absence) dataset using a single draw from the prior distribution, and runs the algorithm on this simulated data; outputs code/output/simulated_data_results.csv
