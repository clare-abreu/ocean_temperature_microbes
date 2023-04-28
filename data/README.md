This folder contains the data used in Abreu & Dal Bello et al (Science Advances, 2023). Information on the datasets can be found the Methods section of the manuscript.

Nicknames used for the datasets in the following files are:

Linnaeus Microbial Observatory: #LMO

ANT 28-5 cruise: #AOT

TARA Oceans project: #TARA

Pivers Island Coastal Observatory (PICO) site: #Ward

Service d’Observation du Laboratoire Arago (SOLA) sampling station: #Med

USC Microbial Observatory at the San Pedro Ocean Time-series station in the San Pedro Channel: #SPOT

P15S GO-SHIP transect: #SPT


The included files are:

# All_data.csv: 
This file contains all estimates of weighted mean copy number (MCN) and weighted mean growth rate (MGR) (obtained by a codon usage bias method), in addition to all metadata available from each dataset. The Estimate column gives the MCN or MGR estimate for each sample.

COLUMNS:

Dataset: official name of dataset

Sample_ID: sample IDs given in datasets; a suffix denoting the dataset nicknames (as described above) is added

Dataset_type: 'time' for time series and 'space' for latitude- and depth-spanning data

Date, Latitude, Longitude, Depth: self-explanatory

Filter_fraction: "FL" for free-living, "SPA" for small-particle-attached, "LPA" for large-particle-attached, and "NF" for no filter

SAR11: indicates whether SAR11 bacteria are included in the estimate

Growth_measure: indicates whether the estimate is of weighted mean copy number (WMCN) or weighted mean growth rate (WMGR)

Copiotrophs: indicates whether bacteria with generation time < 5 hours (as estimated by the codon usage bias method) are included in the estimate

Phototrophs: indicates whether cyanobacteria were included in the estimate

Estimate: an estimate of either WMGR or WMCN

Temperature, Phosphate, Nitrate, Ammonium, Chlorophyll, Nitrite, Salinity, Oxygen, pH, Nitrogen_dioxide, Silicate, Day_length: self-explanatory

DOC/DIC/POC: dissolved organic carbon / dissolved inorganic carbon / particulate organic carbon

Insolation: solar radiation

SAR11_RA: relative abundance of SAR11 bacteria

CN1_RA: relative abundance of bacteria with rRNA copy number = 1

UNITS OF ENVIRONMENTAL METADATA:

Depth: meters

Day length: hours

Temperature: Celsius

Phosphate, Ntirate / Nitrite, Ammonium, DOC (dissolved organic carbon), other nutrients:: micrograms/liter

Salinity: parts per thousand


# Mega_Table.csv:
This file contains all taxa found in all datasets of free-living bacteria, with the exception of the LMO dataset. Since the LMO dataset is so large, it may be accessed on Dryad: https://datadryad.org/stash/share/LLaE8Wv1Wmo6f0M0nlL7S7U2kIfdj1OH9LHfnmKv8mA

COLUMNS:

Phototroph: indicates if the taxa is a phototroph (1) or heterotroph (0) 

Copy Number: the rRNA copy number assigned to each taxa, based on a weighted estimate from the highest-resolution assignment available in the Ribosomal RNA Operon Copy Number Database

Copy Number Classification Level: indicates the taxonomic level of the assignment 

The columns to the right of the taxonomic data are abundance data received directly from published studies, and a suffix/prefix in each column name indicates the dataset. Most datasets use raw counts, but others (i.e. SPOT) use frequencies, but copy number calculations were done each dataset separately, with this difference in mind.


# Mega_Table_metadata.csv:
This file contains some of the data also present in All_data.csv, as well as some WMCN and WMGR calculations for various scenarios, given for each sample.

COLUMNS:

Sample: sample IDs given in datasets; a suffix denoting the dataset nicknames (as described above) is added

WMGR: weighted mean growth rate of the sample, as calculated with the codon usage bias method

WMGR, no SAR11: weighted mean growth rate, not including SAR11 bacteria

WMGR, copio: weighted mean growth rate, for only "copiotrophs," or taxa with generation time < 5 hours

WMGR, Heterotrophs: weighted mean growth rate, not including phototrophs (cyanobacteria)

WMCN: weighted mean copy number of the sample, as calculated with the rrnDB estimates

WMCN, no SAR11: weighted mean copy number, not including SAR11 bacteria

WMCN, no CN1: weighted mean copy number, not including all bacteria with rRNA copy number = 1 (SAR11 is only one example of this)

WMCN, Heterotrophs: weighted mean copy number, not including phototrophs (cyanobacteria)

SAR11 rel. abund.: relative abundance of SAR11 bacteria in the sample

CN1 rel abund: relative abundance of all bacteria with rRNA copy number = 1 in the sample

Temperature, Phosphate, Nitrate / Nitrite / Nitrogen, Chlorophyll, Date, Latitude, Longitude, Depth, Salinity, Oxygen, pH, Ammonia, Nitrogen Dixoide, Day Length, Silicate: self-explanatory

DOC/DIC/POC: dissolved organic carbon / dissolved inorganic carbon / particulate organic carbon

SIOH4: orthosilicic acid

Size fraction lower threshold: The size of the filter, in micrometers, the sample was filtered and only contents above this size were included

Size fraction upper threshold: The size of the filter, in micrometers, the sample was filtered and only contents below this size were included

Filter Fraction: The size of the filter, in micrometers, through which the sample was filtered. Some datasets only include this column, and "3-0.2" means contents between these two sizes were used

Sample Site: A name given by the generators of the data to each site

Insolation: Solar radiation

UNITS OF ENVIRONMENTAL METADATA:
Depth: meters
Day length: hours
Temperature: Celsius
Phosphate, Ntirate / Nitrite, Ammonium, DOC (dissolved organic carbon), other nutrients: micrograms/liter
Salinity: parts per thousand

# rrnDB-5.6_pantaxa_stats_NCBI_WITH_SAR11_CLADES.csv:
This file was downloaded from the Ribosomal RNA Operon Copy Number Database (rrnDB version 5.6, https://rrndb.umms.med.umich.edu/). It contains the rRNA copy number assignments based on taxonomy. We made one edit to this file: the last 12 rows (6821-6832) were added to include SAR11 clades that appeared in datasets and were not represented in the table. These clades are all known to have rRNA copy number = 1.

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
