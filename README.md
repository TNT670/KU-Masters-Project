# KU-Masters-Project
Repository to store the files and code of my M.S. Computer Science project.

## Dependencies/packages required

### R:
Bioconductor (primary resource for most of the microarray analysis tools). Once installed, use BiocManager::install to install:
- affy
- biomaRt
- WGCNA
- preprocessCore
- sva

Also install:
- ggplot2
- ggfortify

### Python
Using a virtual environment either via python or conda can be helpful. A conda virtual environment makes setting up Keras to train MLP models via GPU much easier as well.

Otherwise, be sure the following are installed:
- python3 (3.9+ recommended)
- numpy
- pandas
- matplotlib
- scikit-learn
- tensorflow
- keras
- keras-tuner
- xgboost

Recommended:
- sklearnex (can speed up certain scikit-learn operations if using an Intel machine)
- joblib (to store trained models and skip retraining on re-runs)

## Steps to Run
1. Gather the data. The report uses the following individual datasets:
 - For breast cancer, [GSE5460](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5460), [GSE10780](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10780), [GSE26457](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26457), [GSE29431](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29431), [GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568), and [GSE66162](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66162).
 - For lung cancer, [GSE19188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188), [GSE19804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804), [GSE27262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27262), and [GSE51024](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51024).
 - For colon cancer, [GSE23878](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23878), [GSE71571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71571), [GSE92921](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92921), and [GSE143985](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143985).
 - For simplicity, make a new folder inside the corresponding cancer type folder and name it the GSE ID of the dataset you are downloading. For example, create a "5460" folder inside of the "breast_cancer" folder and place your downloaded CEL files from GSE5460 there, and so on.
2. Copy and paste the R script into the cancer type folder and make any necessary changes as specified in the script.
3. Run the commands in the script. A table will be output at the end; you can specify the name of the table if necessary.
4. Next, open the Jupyter Notebook and ensure that the global variables and the name of the table just output are correct.
5. Observe the order of the datasets and samples in the train dataframe, and adjust the label array so that the corresponding labels are in the same order as the samples if necessary. Changes can be made to any of the models if necessary.
6. Run the Jupyter Notebook cells.
