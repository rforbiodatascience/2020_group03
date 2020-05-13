# Toolbox for analysis and prediction of protein and peptide variant effects
![Image description](https://github.com/rforbiodatascience/2020_group03/blob/master//results/04_heatmaps/heatmap_data_set_score_1.png| width=100))
## 22100 - R for Bio Data Science - Spring 2020
### Group 03 - Begoña Bolós, Jakob Kofoed, Felix Pacheco and Laura Sans
In this toolbox we explored, analyzed and modelled four different deep mutational scanning studies.


## Description

### Introduction
The toolbox found in this repository contains an analysis and predictions tools on four different mutational scanning studies. The aim of this tool box is to predict the activity of a given sequence with or without mutations. The analysis and models proposed, predict the score of the given proteins : BRCA1, ERK2, DLRAP1 and Pab1.


### Workflow
The prediction relies exclusively on the sequence of the protein. To do so, the workflow of the toolbox can be found below.

![](https://github.com/rforbiodatascience/2020_group03/blob/master/doc/external_figures/flowchart.png)<!-- .element height="50%" width="50%" -->

### Repository Structure
The structure of the toolbox can be found below.
![Image description]https://github.com/rforbiodatascience/2020_group03/blob/master/doc/external_figures/00_project_organisation.png| width=100)

### Cleaning and sequence encoding

First each dataset is cleaned to get a "tidy" format where only the valuable data is kept. Then, the sequence strings are edited to reduce the length the protein. The resulting string is encoded to be computed on predictive models such as neural networks. For the encoding of the protein, the use of ten different matrices is implemented. Although, the user is free to use other matrices if the format matches the matrices found in this tool box.


### Predictive models
Once the protein is encoded, three models are available to use; ANN2, ANN (keras) and an elastic net (gml2).

### Shiny App
A shiny app has also been created to perform interactive peptide predictions. The url to the shiny app corresponds to the following link : https://felix-pacheco.shinyapps.io/peptide_score_predictor/

The shiny App only runs for z-scales matrix on the elastic net models for all the proteins. This is due to technical reasons, every prediction relies on the computation of a predictive model. The combination of four datasets, all the encoding options and the three possible models made it difficult to implement all the possible combinations of predictive models.

Example sequence for Pab1 prediction: "ANLHPDIDNKALYDTFSVFGDLLSSKIATDENGKSKGFGFVHFEEEGAAKEAIDALNGMLLNGQEIY"
It must be noted that, for the prediction, the string has to be the same length as the sequence used to train the model.

## Installation



## Usage 

The analysis tool can be used by running the 00_doit.R file. The running can be customized by the user, e.g to only run one model.

## License

MIT





## Data sources:

#### Data set 1 : BRCA1
* reference: https://doi.org/10.1534/genetics.115.175802


#### Data set 2 : ERK2
* reference: https://doi.org/10.1016/j.celrep.2016.09.061


#### Data set 3 : DLRAP1
* reference: https://www.mavedb.org/scoreset/urn:mavedb:00000036-a-1/


#### Data set 4 : Pab1.
* reference: https://doi.org/10.1371/journal.pgen.1004918

