# ml-sleep
Script used to build machine learning model to predict sleep genes. Features are built using two lines of information. <br />
 -  Genome wide data with evidence of sleep genes enrichment.
 -  Similarity index of this gene to sleep genes given the molecular context.

Example inputs needed to run the ML model are prepared in DATA. 
 -  Step 1. Download human and mouse gene info from https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/, and put these two files in REFERENCE. <br />
 -  Step 2. Run "./R/process_rawInput.R" to prepare input table. <br />
 -  Step 3. Run "./R/calc_EvidenceFactors.R" to calculate evidence factors for each genome-wide data. New table will only keep data with high evidence factors. <br />
 -  Step 4. Run random forest prediction using RF_prediction.py <br />


