This folder contains files related to a test run of Badsha's code for 5 tissues after changing the pc selection format

*discovered that badsha's code preformed a pvalue adjustment by the default method in corr.test() before being passed to 
 qvalues() leading to an extremely conservative selection of PC's for each trio

we change the 'adjust' parameter of corr.test() to 'adjust="none" and then carry out the PC selection

the five tissues are:

WholeBlood
AdiposeSubcutaneous
ArteryTibial
MuscleSkeletal
SkinSunExposed

the files 

./List.significant.asso1.tissue.RData

contain the list of selected PC's for each trio under this adjustment to the code in the form:
             trio1     trio2   ....     trioN
[[1]] PC1     


./List.Match.significant.trios.tissue

contains the list in which the list elements are the trios and the vectors are integer column indexes of the PC matrix