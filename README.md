# KU-Masters-Project
Repository to store the files and code of my M.S. Computer Science project.

## Steps to Run
1. Gather the data. The report uses the following individual datasets:
 - For breast cancer, [GSE5460](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5460), [GSE10780](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10780), [GSE26457](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26457), [GSE29431](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29431), [GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568), and [GSE66162](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26457).
 - For lung cancer, [GSE19188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188), [GSE19804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804), [GSE27262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27262), and [GSE51024](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51024).
 - For colon cancer, [GSE23878](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23878), [GSE71571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71571), [GSE92921](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92921), and [GSE143985](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143985).
 - For simplicity, make a new folder inside the corresponding cancer type folder and name it the GSE ID of the dataset you are downloading. For example, create a "5460" folder inside of the "breast_cancer" folder and place your downloaded CEL files from GSE5460 there, and so on.
2. Copy and paste the R script into the cancer type folder and make any necessary changes as specified in the script.
3. Run the script. A table will be output.
