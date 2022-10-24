# Final project for the "Python for biologists" course at Gothenburg University

This is a script for a basic statistical description of data about gene expression levels.

## Description

The script can be divided into 6 parts:
1. Standardization of values
2. Principal component analysis (PCA) 
3. T-test of every marker between two groups
4. P-value adjustment using 3 different ways
5. Volcano plot
6. Heatmap

The goal of the script is two explore the differences between the two groups of people in the dataset and to detect the markers driving these differences.

## Installation

At the top of the script are mentioned the libraries need for running the script (Pandas, Sklearn, Matplotlib, Seaborn, bioinfokit, statsmodels, matplotlib-venn, adjusttext, tabulate, textwrap3). Please make sure you have installed all of them, otherwise you can install them in your conda environment with:

   conda install pandas
   conda install -c intel scikit-learn
   conda install -c anaconda statsmodels 
   conda install -c bioconda bioinfokit 
   conda install -c conda-forge matplotlib-venn
   conda install -c conda-forge adjusttext
   conda install tabulate
   conda install -c conda-forge textwrap3
   conda upgrade bioinfokit

## Usage

The script was created in Spyder 5.1.5 using Python 3.9.
Atfer having cloned or downloaded the data on your local computer, the script can be run in an IDE like Spyder or in the Anaconda
Prompt with the command:

    python genes_expression_project.py. 
    
This will generate 3 plots in your working folder: one for the PCA, one volcano plot and one heatmap. 

## Contributions

I would like to thank the teachers and the other course participants for helping me with this project with their comments and assisting me with the script.

