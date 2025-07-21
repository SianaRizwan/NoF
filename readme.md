
<!-- 1. For 0-1 knapsack (profit-maximization):
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


The above process should be enough to generate the fuinal results. -->
## 📘 0-1 Knapsack ML Pipeline Guide

This document explains the complete process to conduct experiments on the 0-1 Knapsack Problem, for both profit-maximization and profit-minimization variants. 

🛠️ All scripts were tested on Linux with Python 3.12+.
## 📂 Folder Structure Overview
```
NoF/
├── Analysis/                                  # Analysis scripts and metrics
│   ├── llm_eval/
│   └── top_feature_per_metric/
├── Codes/                                     # Codebase for the pipeline
│   ├── Algorithms/
│   │   ├── profit_maximization/               # Algorithms/Solvers for maximization
│   │   └── profit_minimization/               # Algorithms/Solvers for minimization
│   ├── Feature/
│   │   ├── extract_features_kp.cpp            # Feature generator (C++)
│   │   └── optimality_gap.ipynb               # Optimality gap vs Gurobi
│   ├── Machine Learning Models/
│   │   ├── Classification/                    # ML models for classification
│   │   └── Regression/                        # ML models for regression
│   └── Others/
│       ├── extract_bins.ipynb                 # Binning instances by performance metric
│       ├── merge_algo_feature.ipynb           # Merges solver outputs with extracted features
│       ├── merge_files.ipynb                  # Combines multiple solver runs into a single CSV per algorithm
│       └── train_test_val_split.ipynb         # Creates train-test-validation splits for ML model training
├── Dataset/                                   # Input problem instances
│   ├── maximization/
│   │   ├── knapsack_dataset.csv
│   │   └── training_data/                     # Processed data (bins, splits)
│   └── minimization/
│       ├── knapsack_dataset.csv
│       └── training_data/                     # Processed data (bins, splits)
├── Results/                                   # Model outputs and results
│   ├── Maximization/
│   │   ├── Classification/
│   │   ├── Ensembles/
│   │   └── Regression/
│   └── Minimization/
│       ├── Classification/
│       └── Ensembles/
└── Framework.pdf                              # High-level framework overview
```

## ✅ Steps to Run the Experiment
### Step 1: Solve Knapsack Instances:
- Go to Codes/Algorithm/profit_maximization/ .
- Run each algorithm on the dataset file knapsack_dataset.csv located in Dataset/maximization/ folder.
- Set RAM and CPU manually inside the code.
- Save outputs in separate folders for each algorithm.

### Step 2: Feature Extraction
- Compile and run extract_features_kp.cpp (inside Codes/Feature/) to generate instance-wise features.

### Step 3: Merge Feature and Algorithm Output
- Run merge_algo_feature.cpp (inside Codes/Others/) to merge algorithm outputs with features.
- Run merge_files.ipynb (inside Codes/Others/) to combine results from different runs into one final CSV file for each algorithm.

### Step 4: Calculate Optimality Gap
- Run optimality_gap.ipynb (in Codes/Feature/) to compute optimality gap against Gurobi for each instance.

### Step 5: Create Bins for ML Classification
- Run extract_bins.ipynb (inside Codes/Others/) specifying path to the final merged CSV for each algorithm to generate bins for performance metrics. 

### Step 6: Split Dataset into Train/Test/Val
- Run train_test_val_split.ipynb (in Codes/Others/). This will generate three files for each algorithm: _train.csv, _test.csv, _val.csv.

### Step 7: Train ML Models
- Go to the folder Codes/Machine Learning Models.
- For running each of the classification-based scripts, folder paths containing bin files and train-test-val files need to be provided (e.g., path to bin files: Dataset/maximization/training_data/bins and path to train-test-val files: Dataset/maximization/training_data/algorithms)
- For running regression-based scripts, folder path containing the train-test-val files need to be provided only (e.g., path to train-test-val files: Dataset/maximization/training_data/algorithms)
- Run the ML scripts to train and generate predictions.

📂 Example folder structure for training data is provided in: Dataset/maximization/training_data/ 

## 📌 Notes on the 0-1 knapsack Variants

### Profit-Maximization
- Dataset path: Dataset/maximization/knapsack_dataset.csv
- Solver codes: Algorithm/profit_maximization/

### Profit-Minimization
- Dataset path: Dataset/minimization/knapsack_dataset.csv
- Solver codes: Algorithm/profit_minimization/
- All other steps remain exactly the same.