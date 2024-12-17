The files in this folder are as follows:

files of the form:

'tissue.name'.eGenes.V8.unique.RData

Are intermediated files that are needed for use in the functions contained in ADDIS_rerun.R
These files are created from the function get.unique.genesfile() in ADDIS_rerun.R
The files were saved to this folder for to eliminate their compilation step in rerun.ADDIS()

files of the form:

List.models.'tissue.name'.all.RData

Are 5 element lists for each of the 5 model types giving the index of the trios classified to each model

i.e.

M0 --> 1, 5, 7, 19, ...

M1 --> 2, 4, 9, 20, ...


These .RData lists are used in the functions contained in the R-scripts found in the folder(s) ./ADDIS_verify and ./GMACanalysis


files of the form

All.models.'tissue.name'.all.RData

are summarizing arrays showing the total number of trios classified to each of the model types (including the type 'other') 


The file ./tissuenames.csv

contains the 3 tissue name formats commonly associated with both GTEx and Badsha's files. It is used in several of the subdirectories: dim 48 X 3

The file ./All.models.result.csv

is a summarizing table combining the results from all of the All.models.'tissue.name'.all.RData files for all tissues

the files 

./List.significant.asso1.tissue.RData

contain the list of selected PC's for each trio under this adjustment to the code in the form:
             trio1     trio2   ....     trioN
[[1]] PC1     


./List.Match.significant.trios.tissue

contains the list in which the list elements are the trios and the vectors are integer column indexes of the PC matrix