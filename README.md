# ml-sleep
Script used to build machine learning model to predict sleep genes. Features are built using two lines of information. <br />
  1. Genome-wide datasets with evidence of sleep genes enrichment. <br />
  2. Gene and protein knowledge from annotated gene set collections are scored using Jaccard index (JI). <br />

<br />

**Example inputs to run the ML model are prepared in 'DATA'. Follow the steps to prepare and process data as needed.**
 -  Step 1. Download human and mouse gene info from https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/. Move these two files into 'REFERENCE'. <br />
 -  Step 2. Run the code in "./R/process_rawInput.R" to prepare input table. <br />
 -  Step 3. Run the code in "./R/calc_EvidenceFactors.R" to calculate evidence factors for each genome-wide data. Five processed datasets are available as example. The new table will only remove those data with evidence factors less than 3. <br />
 -  Step 4. Run random forest prediction using "./pycharm/RF_prediction.py". <br />
 -  Step 5. Run "./R/process_RF.prediction.R" to summarize predictions made by the ML model. <br />


