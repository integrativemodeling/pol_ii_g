# pol_ii_g
These scripts demonstrate the use of IMP, MODELLER, and PMI in the modeling of the POl II (G) complex using DSS chemical cross-links.

Localizing gdown1 on POL II using integrative structure modeling

Input data (directory data):
40 DSS cross-links involving gdown1 were identified via mass spectrometry, 14 of these cross-links were intramolecular and 26 were intermolecular with Pol II. The Pol II structure used was obtained from the PDB (code 5FLM); it was determined primarily based on a cryo-EM density map at 3.4 Å resolution (EMDB code: 3218) {Bernesky, 2016 Nature}.
Representation of gdown1 relied on (i) secondary structure and disordered regions predicted by PSIPRED based on the gdown1 sequence (Buchan et al., 2013; Jones, 1999).

Running the simulations (directory production_scripts)
-sample.py: pmi modeling scripts for running the production simulations: The search for good-scoring models relied on Gibbs sampling, based on the Metropolis Monte Carlo algorithm. We suggest producing at least 5 million models from 100 independent runs, each starting from a different initial conformation of gdown1 to have proper statistics.
-submit.sub: SGE-cluster based submission scripts to run automatically 100 independent runs.
-The compressed 100 independent trajectories are accessible at: X.tar

Analyzing the simulations (directory analysis_scripts)
various scripts to analysis the simulations. we give more details scripts that allows us to test for sampling convergence.
-select_good_scoring_models.py: python script that read and select RMF files based on good scoring model criteria explained in the methods section.
-random_subsets_best_score_convergence.py: we test if adding more models improve our sampling of top scores. input file to the script is a list of all the Total_Score from the simulations. For each subset, we perform 100 sub-samplings to compute error bars. We give the file twice to check how the error bars varies.
-Kolmogorov_Smirnov_2Samples.py: we test if the distribution of from two independent subsets are not unsimilar.

-The good-scoring models that have been selected for precision-based clustering based on RMSD metric are located in: XXX 
The file name format is ${trajectory_number}_${frame_number}.rmf3

-sample_precision.py: Given a RMSD matrix, we compute the 𝛘2-test for homogeneity of proportions{McDonald, 2014} and compute the best precision for which sampling has converged.
-rmsf.py: Given a list of structures, compute the average RMSF, which indicates the precision of the structures in a cluster.

Plotting the results (directory plotting_scripts)
Python scripts for generating figures from the paper.


Publication: 
Architecture of Pol II(G) and molecular mechanism of transcription regulation by Gdown1 

Miki Jishage, Xiaodi Yu, Yi Shi, Sai J. Ganesan, Andrej Sali, Brian T. Chait, Francisco Asturias, and Robert G. Roeder
https://doi.org/10.1038/s41594-018-0118-5
