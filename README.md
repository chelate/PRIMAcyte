# AkerselvaCytometry
A collection of Julia notebooks designed to get your IMC analysis off and running.
Assumes that `.mcd` files  have already been unpacked by the Steinbeck docker tool and segmented using deepcell/mesmer or CellPose.

This repository contains a set of Julia Pluto notebooks that can be run in succession to perform data analysis.

## How to Run

1. Clone the repository:
	In terminal
	```
	cd pathto/SteinbockProject/
	git clone https://github.com/chelate/AkerselvaCytometry
	cd AkerselvaCytometry
	```
	
2. Install Julia and Pluto.jl (if not installed):
	- Julia: follow instructions at	
		https://julialang.org/downloads/
	
	- Pluto: open Julia in terminal and type
	```
	] add Pluto
	```
	
3. Run the analysis pipeline:
	from pathto/SteinbockProject/AkerselvaCytometry/
	- To run them all at once
	```
	julia run_analysis.jl
	```
	
	- To run them notebook by notebook and generate images, in julia terminal at `pathto/SteinbockProject/AkerselvaCytometry/`
	julia```
	using pluto
	Pluto.run()	
	```
	
## Notebooks
- `01_QualityControl.jl`: Cleans the raw data, removes low quality cells and images, attatches patient identifiers / metadata such as tme block
	Optional: folder called lastik_labels can be used to provide probabilities to filter out low quality regions.
	
- `02_BatchCorrections.jl`: Extracts features from the cleaned data. Generates a new c1.csv in AkerselvaPipeline/data that is appropriate for use in CytoBank or similar.

- `03_CelltypeAnnotations.jl`: Clusters cells and annotates them according to a binary tree. The clustering is saved in AkerselvaPipeline/data for stability.  

- `04_CelltypeAnnotations.jl`: Clusters cells and annotates them according to a binary tree. The clustering is saved in AkerselvaPipeline/data for stability.  