# Dependencies
This program is required to use pymol, Biopython.PDB and other normal packages.
It is written in python3.

# Program components
Our program is used to solve 2 tasks. Task 1 is to analyze relationships within one ensemble , task 2 is to analyze relationships between different ensembles.
The program consists of 5 simple .py file. 

## Task 1
To run the program, first we need to execute "task1.py" file and input the id of ensemble we want to study. Then the pdb files of the ensembles will be downloaded automatically. 
+ The features of one conformation(output a) will be save in the path "/data/output/{PED\_id}/..." .
The feature files are saved as name "feature\_conf\_XXX.json". Here, XXX can be "rg", "asa","ss" and "dm". 
 * "rg" means radius of gyration.
 * "asa" means relative accessible surface area.
 * "ss" means secondary structure.
 * "dm" means distance matrix.
+ The graphs reflecting the distance between conformations(output b) are saved in the path "/data/output/{PED\_id}/..." as name "task1\_output\_b\_{PED\_id}_X.png", while X means different angle of the 3D picture.
+ Then we execute "task1\_output\_c.py" and get the Pymol image of the conformations(output c) we studied. The image is also saved in the path "/data/output/{PED_id}/..." . It is saved as name "pymol\_image.png".

## Task 2
We execute "task2.py" and get the 3 outputs. All of them are saved in the same path before.
+ The features of one ensemble(output a) are saved as name "feature\_ens\_sd\_XXX.json". Here, XXX can be "rg", "ss\_entropy", "med\_asa", "med\_rmsd", "med\_dm" and "sd\_dm".
  * "rg" means radius of gyration for each conformation in the ensemble.
  * "ss\_entropy" means secondary structure entropy for each position across ensemble conformations.
  * "med\_asa" means median solvent accessibility for each position across ensemble conformations.
  * "med\_rmsd" means median RMSD for each position across ensemble conformations.
  * "med\_dm" means median distance of each pair of equivalent positions across ensemble conformations.
  * "sd\_dm" means standard deviation of the distance of each pair of equivalent positions across ensemble conformations.
+ The dendrogram and heatmap(output b) are saved as name "dendrogram\_distance\_{PED\_id}.png" "heatmap\_distance\_{PED\_id}.png".
+ The plot showing features values along sequence positions(output c) is saved as name "local\_score\_{PED\_id}.png".
