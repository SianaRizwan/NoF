## Steps to conduct the experiment ##
1. For 0-1 knapsack (profit-maximization):
- Under folder 'Algorithm/profit_maximization', the algorithm/solver codes are present that need to be run.
- Under folder 'Dataset/maximization', the knapsack_dataset.csv contains the instances that is to be solved by the algorithms.
- Go to 'Algorithm/profit_maximization' to run each algorithm on 'Dataset/maximization/knapsack_dataset.csv'. You need to specify the RAM and CPU in the codebase, manually.
- Save outputs in separate folders for each algorithm.
- Run extract_features_kp.cpp in 'Codes/Feature' to generate instance features.
- Run merge_algo_feature.cpp in 'Codes/Others' to combine features with algorithm outputs.
- Run merge_files.ipynb in 'Codes/Others' to merge all outputs per algorithm.
- Run optimality_gap.ipynb in 'Codes/Feature' to calculate the optimality gap for each instance compared to gurobi. At the end there should be one file per algorithm.
- To extract bins for each performance metric per algorithm, run extract_bins.ipynb specifying the path to the final merged csv file for each algorithm.
- To divide train-test-validation dataset, run train_test_val_split.ipynb in 'Codes/Others' folders. This will generate 3 files (_train, _test, _val . csvs) per algorithm.
- The demo of the folder structure required can be found in 'Dataset/maximization/training_data' folder.
- Lastly, to run the ML models, got to the 'Machine Learning Models'. Two folder paths need to be specified; one containing the bins per algorithm ('Dataset/maximization/training_data/bins)' and another containing the train-tes-val files per algorithm ('Dataset/maximization/training_data/algorithms').


2. For 0-1 knapsack (profit-minimization):
- To run the algorithms for this variant, go to folder 'Algorithm/profit_minimization'. And, for dataset, refer to knapsack_dataset.csv under the folder 'Dataset/minimization'.
- The rest of the steps for preparing the dataset and training the models are same as the profit-maximization variant.


The above process should be enough to generate the fuinal results.