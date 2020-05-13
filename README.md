# Toolbox for analysis and prediction of protein and peptide variant effects
![Image description](https://github.com/rforbiodatascience/2020_group03/blob/master//results/04_heatmaps/heatmap_data_set_score_1.png)
## 22100 - R for Bio Data Science - Spring 2020
### Group 03 - Begoña Bolós, Jakob Kofoed, Felix Pacheco and Laura Sans
In this toolbox we explored, analyzed and modelled four different deep mutational scanning studies.


## Description

### Introduction
The toolbox found in this repository contains an analysis and predictions tools on four different mutational scanning studies. The aim of this tool box is to predict the activity of a given sequence with or without mutations. The analysis and models proposed, predict the score of the given proteins : BRCA1, ERK2, DLRAP1 and Pab1.


### Workflow
The prediction relies exclusively on the sequence of the protein. To do so, the workflow followed as well as the structure of the repository can be found below.

![Image description](https://github.com/rforbiodatascience/2020_group03/blob/master/doc/external_figures/flowchart.png)[width=20%]

First each dataset is cleaned to get a "tidy" format where only the valuable data is kept. Then, the sequence strings are encoded to be computed on predictive models such as neural networks. For the encoding of the protein, the use of ten different matrices is implemented. Although, the user is free to use other matrices if the format matches the matrices found in this tool box.

Once the protein is encoded, three models are available to use; ANN2, ANN (keras) and an elastic net (gml2).

### Shiny App
A shiny app has also been created to perform interactive peptide predictions. The url to the shiny app corresponds to the following link : https://felix-pacheco.shinyapps.io/peptide_score_predictor/

The shiny App only supports the use of the elastic net as a predictor, some string examples to run on each protein can be found in the shiny app code as a comment (in the server section).

## Installation


## Usage 


Describe usage of the analysis scripts


## Contributing



## License

MIT





## Data sources:

### Data set 1
* data_set_1 ref: https://doi.org/10.1534/genetics.115.175802
* data_set_1 <- read_excel('./data/_raw/genetics.115.175802-6.xls')

### Data set 1
* data_set_2 ref: https://doi.org/10.1016/j.celrep.2016.09.061
* data_set_2 <- read_excel('./data/_raw/1-s2.0-S2211124716313171-mmc2.xlsx', sheet = 'Supplemental_Table_1')

### Data set 3
* data_set_3 ref: https://www.mavedb.org/scoreset/urn:mavedb:00000036-a-1/
* data_set_3 <- read_csv('./data/_raw/urn_mavedb_00000036-a-1_scores.csv', skip = 4)

### Data set 4
* data_set_4 ref: https://doi.org/10.1371/journal.pgen.1004918
* data_set_4 <- read_tsv('./data/_raw/S1_Table.txt')
