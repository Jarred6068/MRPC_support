Trio Alignment

- to form trios, the UCSC refernce gene name of both expression probes and methylation probes were used. For genes with multiple methylation probes
- for the same gene, each was paired with the respective expression probe to form a trio 

- the final list of methylation probe ID's aligned with expression probes is given in the Rdata file 
		/mnt/ceph/jarredk/Methyl/Used.Mprobes.final.Rdata



-CalcCor meta data files:

- Expression meta data was obtained via BioMart in BSGS to acquire the expression probe coordinate, ILMN_ID, and Chromosome number. The Expression
  matrix was aligned with the biomart data via the ILMN_ID number
		BM file : "/mnt/ceph/jarredk/Methyl/mart_export_GB37.txt"   dim   (260931 X 5)  
		expression file:  "/mnt/ceph/jarredk/Methyl/ExpressData/Peerdata.bsgs.RData"   dim   (859 X 24317)

- The ~~3000 probes that did not match up with biomart are listed by their gene name in the file:

		/mnt/ceph/jarredk/Methyl/unmatched_genes_biomart.txt

		the final files /Expression_data_BM_aligned.Rdata and /meta_E.Rdata in the folder /Correlation_Calculation
		are of the dimension (606 X 20578) and (20578 X 5) respectively
#NOTE:::
	For row alignment of the methylation and expression matrices please see the normalization README file

- Methylation meta data was obtained by using the ID's of probes included in the final list for trio formation and aligning these probes to the 
  to the file /mnt/ceph/megheib/M_G_data/GPL13534_M.csv or /mnt/ceph/jarredk/Methyl/GPL13534_HumanMethylation450_15017482_v.1.1.csv in my directory.
  to extract the 

- The final list of methylation probes used in trio construction is given in the file 
		/mnt/ceph/jarredk/Methyl/Used.Mprobes.final.Rdata

	and contains the ID's for the 363,516 retained probes


