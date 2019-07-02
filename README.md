# CPPNet: Coding Protein Predication Based on Deep Network of RNA Sequences  
A deep learning-based method to predict protein-coding potential of RNA sequences
# Create virtual environments and install dependencies
conda create -n CPPNetEnv python=3.6  
conda activate CPPNetEnv (or source activate CPPNetEnv)  
conda install r  
install R package "LncFinder"  
      install.packages("LncFinder")  
pip install numpy  
pip install pandas  
pip install sklearn  
pip install biopython  
pip install tensorflow  
pip install keras  
pip install rpy2==3.0.1  
# Inputs
Fasta file
# Usage
## Human_CPPNet model
$ ./code/Human_CPPNet.py -i input.fasta -o output.csv
## Cross-species CPPNet model
$ ./code/Cross_species_CPPNet.py -i input.fasta -o output.csv
