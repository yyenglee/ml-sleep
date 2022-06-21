# ml-sleep
Updated pipeline (06/21/22), rewrite the code in Python to increase efficiency.

Script used to build machine learning model to predict sleep genes. Features are built using two lines of information. <br />
  1. Genome-wide datasets with evidence of sleep genes enrichment. <br />
  2. Gene and protein knowledge from annotated gene set collections are scored using Jaccard index (JI). <br />

<br />

**Example inputs to run the ML model are prepared in 'DATA'. Follow the steps to prepare and process data as needed.**
 -  Step 1. Download human and mouse gene info from https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/. Move these two files into 'REFERENCE'. <br />
 -  Step 2. Run the code in example_run_script.py. <br />
            This pipeline includes - <br />
            1. Calculate evidence factors based on the seed genes for each datasets in the genome-wide table. Select those that can provide information (maxEF>3). <br />
            2. Additional features based on gene sets annotation information and Jaccard index. <br />
            3. To decide if going to run the <br />
                (i) model evaluation with seed genes, <br />
                (ii) model evaluation with randomly assign labels, <br />
                (iii) predictions using random forest models. <br />
              Whether to run each of these steps is defined at the beginning of the pipeline. <br />
