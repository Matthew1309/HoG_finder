This repository contains my work in tetramer nucleotide frequency dimension reduction. The purpose of these scripts is to take any given fasta input, chop it up, encode each fragment as a 256 dimension vector, then graph/cluster them.

Currently it is 8/28/2022, and the last time I worked on this was in 10/2021, and I didn't make a readme or conda env then.

I believe the workflow is to use ```genome_to_tetramer.ipynb``` to turn a given .fa file into a csv table where rows are fragements and columns are tetramers. I have some clustering tutorials that I followed, but I want to start fresh now that I've used Seurat.

Input: csv of tetranucleotide frequencies
Step1: Normalize the columns into Z-scores
Step2: Perform PCA
Step3: Peform t-SNE on PCs rather than raw data
Step4: Find clusters using any method I want! K-means, heiarchical, or whatever seurat does.

I will put this into a file called `HoG_workflow.ipynb`
I also added the ```stats.yml``` which should set up the proper conda env with the command `conda env create -f "stats.yml"`.