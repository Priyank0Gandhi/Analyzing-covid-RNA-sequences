# Analyzing-covid-RNA-sequences

I. Introduction:
--
In this project, I have retrieved and analyzed the RNA sequence data of COVID's prominent strains, Delta and Omicron. 
RNA, a nucleic acid found within every living cell, comprises a single strand composed of diverse arrangements of four nucleotides: uracil, cytosine, adenine, and guanine. 
This RNA serves as the fundamental blueprint for COVID, facilitating cellular entry and self-replication of the virus.

II. Dataset:
--
Using the available data from NIH (National Institutes of Health), the federal government agency in the U.S.
Documentation of the data is available here: https://www.ncbi.nlm.nih.gov/datasets/docs/v1/data-packages/sars-cov-2-genome/
The drive link for the dataset can be found here: https://drive.google.com/file/d/1S2ZDjdRkY78kZxBtc9YNUh0mByTHXQ23/view

The Dataset consisted of 847,791 records and  17 attributes.

III. Data Preparation:
--
Firstly, the metadata of Covid RNA Sequences was loaded in the Pandas data frame.
Then, the datatypes of each column in the dataset were checked for appropriateness. "Collection Date" column into the appropriate date datatype. And the datatypes of the rest of the columns seemed suitable for analysis.
Further, the "Geolocation" column was split to obtain the continent name and country name in two separate columns.

IV. Exploratory Data Analysis:
--
1.) At what point in time was the first COVID RNA sequence gathered from each of the continents?

--> 
![image](https://github.com/Priyank0Gandhi/Analyzing-covid-RNA-sequences/assets/96395339/e9c17e60-8bc1-42df-8d89-d8cdcb43112c)

2.) How many sequences were collected on each continent?

--> 
![image](https://github.com/Priyank0Gandhi/Analyzing-covid-RNA-sequences/assets/96395339/28982409-ff2a-4e12-a928-19150e5e4600)            ![image](https://github.com/Priyank0Gandhi/Analyzing-covid-RNA-sequences/assets/96395339/9eb2d24a-7f0f-40e0-ac5e-395801b0e1e3)

3.) How long are the shortest and longest sequences? Look at the outliers, if any, to estimate their representation.

--> 
![image](https://github.com/Priyank0Gandhi/Analyzing-covid-RNA-sequences/assets/96395339/ec57a92a-d69a-43df-aebd-8963f8ca8fd5)
    Values less than 29639 can be considered outliers.
  Filtered out all the outlier records from the data frame.
    
4.) How many samples were collected by month? Are there any trends?

--> 
![image](https://github.com/Priyank0Gandhi/Analyzing-covid-RNA-sequences/assets/96395339/5c44efe1-4f24-4fde-9a2a-bae0ececf745)

 V. Identifying the RNA sequences to be downloaded, for analysis:
--







