This folder contains datasets formatted to be input into the functions in calc_WMCN_WMGR_functions.py. Information on these datasets is available in the manuscript's Methods section. The only datasets not available here are LMO and SPT, which are very large and available on Dryad (https://datadryad.org/stash/share/LLaE8Wv1Wmo6f0M0nlL7S7U2kIfdj1OH9LHfnmKv8mA ).

Nicknames used for the datasets in these files are:

Linnaeus Microbial Observatory: LMO

ANT 28-5 cruise: AOT

TARA Oceans project: TARA

Pivers Island Coastal Observatory (PICO) site: Ward

Service d’Observation du Laboratoire Arago (SOLA) sampling station: Med

USC Microbial Observatory at the San Pedro Ocean Time-series station in the San Pedro Channel: SPOT

P15S GO-SHIP transect: SPT


LMO and AOT datasets have multiple filter fractions; the files containing "free" in their names include only free-living filter fractions.

COLUMNS:

OTU_ID: a generic identifier; has no meaning across datasets

Phylum, Class, Order, Family, Genus, Species, Sequence: when available for each taxon

Phototroph: Indicates whether the taxon is phototrophic (1) or heterotrophic (0)

Copy Number: The estimated rRNA copy number for this taxon

Copy Number Classification Level: Indicates which taxonomic classification level was used for assigning rRNA copy number (the highest-resolution level available in the rrnDB)

Other columns: All samples from the given dataset. Numbers in these columns represent sequenced abundances.

