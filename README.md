# Toolbox for analysis and prediction of protein and peptide variant effects
![Image description](https://github.com/rforbiodatascience/2020_group03/blob/master/doc/front.png)
## 22100 - R for Bio Data Science - Spring 2020
### Group 03
In this toolbox we explored, analyzed and modelled four different deep mutational scanning studies.

### Data sources:

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
