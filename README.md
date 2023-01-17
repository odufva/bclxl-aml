# bclxl-aml
Scripts to reproduce figures and analyses in the manuscript "Erythroid/megakaryocytic differentiation confers BCL-XL dependency and venetoclax resistance in acute myeloid leukemia" Kuusanmäki, Dufva et al. Blood 2022

## To reproduce the results, obtain source data from Synapse:
- Get synapse credentials https://www.synapse.org
- Access synapse project syn24200411 (https://www.synapse.org/bclxl_aml)
- Download project data:
	- Input files individually (see scripts for filenames and download from SYNAPSE) (Recommended) 
	- Programmatic access (synapse, check synapse IDs from synapse):
		```
		pip install synapseclient
		synapse get synapseID
		```
	- Synapse bulk (SIZE):
		```
		pip install synapseclient
		synapse get syn24200411 -r
		```
		
		
## Processed scRNA-seq data is also available at ArrayExpress
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12607
