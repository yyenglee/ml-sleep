# ml-sleep
Script used to build machine learning model to predict sleep genes. Features are built using two lines of information. <br />
 -  Download human and mouse gene info from https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/, and put these two files in REFERENCE. <br />
 -  Run "./R/process_rawInput.R" to prepare input table. <br />
 -  Run "./R/calc_EvidenceFactors.R" to calculate evidence factors for each genome-wide data. New table will only keep data with high evidence factors. <br />
 -  Run random forest prediction using RF_prediction.py <br />


