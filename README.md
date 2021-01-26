# bclxl-aml
Scripts to reproduce figures and analyses in the manuscript "Erythroid/megakaryocytic differentiation confers BCL-XL dependency and venetoclax resistance in acute myeloid leukemia" Kuusanmäki, Dufva et al. bioRxiv 2021

## To reproduce the results, obtain source data from Synapse:
- Get synapse credentials https://www.synapse.org
- Access synapse project SYNAPSE
- Download project data:
	- Input files individually (see scripts for filenames and download from SYNAPSE) (Recommended) 
	- Programmatic access (synapse, check synapseID_Filename.txt for accession codes):
		```
		pip install synapseclient
		synapse get synapseID
		```
	- Synapse bulk (SIZE):
		```
		pip install synapseclient
		synapse get SYNAPSE -r
		```
