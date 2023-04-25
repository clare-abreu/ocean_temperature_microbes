This folder contains the data used in Abreu & Dal Bello et al (Science Advances, 2023).

The included files are:

# All_data.csv: 
This file contains all estimates of weighted mean copy number (MCN) and weighted mean growth rate (MGR), in addition to all metadata available from each dataset. The Estimate column gives the MCN or MGR estimate for each sample, with the Growth_measure column indicating which type of estimate this number is. The Filter_fraction column indicates "FL" for free-living, "SPA" for small-particle-attached, "LPA" for large-particle-attached, and "NF" for no filter. The Copiotrophs and Phototrophs columns indicate whether copiotrophs (rRNA copy number > 1) and phototrophs (cyanobacteria) were included in the estimate. The SAR11-RA and CN1_RA columns indicate the relative abundance of SAR11 and all bacteria with rRNA copy = 1 in each sample.


# Mega_Table.csv:
This file contains all taxa found in all datasets of free-living bacteria, with the exception of the LMO dataset (which is by far the largest dataset, and may be accessed on Dryad: https://datadryad.org/stash/share/LLaE8Wv1Wmo6f0M0nlL7S7U2kIfdj1OH9LHfnmKv8mA ). The Phototroph column indicates if the taxa is a phototroph (1) or heterotroph (0). The Copy Number column indicates the rRNA copy number assigned to each taxa, based on a weighted estimate from the highest-resolution assignment available in the Ribosomal RNA Operon Copy Number Database (rrnDB version 5.6, https://rrndb.umms.med.umich.edu/), and the Copy Number Classification Level column indicates the taxonomic level of the assignment. The columns to the right of the taxonomic data are abundance data received directly from published studies, and a suffix/prefix in each column name indicates the dataset. Most datasets use raw counts, but others (i.e. SPOT) use frequencies, but copy number calculations were done each dataset separately, with this difference in mind.


# Mega_Table_metadata.csv:
This file contains some of the data also present in All_data.csv, as well as some WMCN and WMGR calculations for various scenarios, given for each sample.


# rrnDB-5.6_pantaxa_stats_NCBI_WITH_SAR11_CLADES.csv:
This file was downloaded from rrnDB. It contains the rRNA copy number assignments based on taxonomy. We made one edit to this file: the last 12 rows (6821-6832) were added to include SAR11 clades that appeared in datasets and were not represented in the table. These clades are all known to have rRNA copy number = 1.


# LMO Accession numbers:
The LMO 16S data accession numbers are:

PRJEB52828          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB52828

PRJEB52855          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB52855

SRP048666            https://www.ebi.ac.uk/ena/browser/text-search?query=SRP048666

PRJEB52837          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB52837

PRJNA260662        https://www.ebi.ac.uk/ena/browser/view/PRJNA260662

PRJEB42455          https://www.ebi.ac.uk/ena/browser/view/PRJEB42455

PRJEB52851          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB52851

PRJEB52854          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB52854

PRJEB52627          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB52627

PRJEB52772          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB52772

PRJEB52496          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB52496

PRJEB52780          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB52780

PRJEB52782          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB52782

PRJEB52850          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB52850

PRJEB56744          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB56744

PRJEB56745          https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB56745
