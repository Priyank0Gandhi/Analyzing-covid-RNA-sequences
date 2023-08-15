# Analyzing-covid-RNA-sequences

## Table of Contents:
I.[Introduction](#Introduction)  
II.[Dataset](#Dataset)  
III.[Data_Preparation](#Data_Preparation)  
IV.[Exploratory_Data_Analysis](#Exploratory_Data_Analysis)  
V. [Selecting_RNA_for_analysis](#Selecting_RNA_for_analysis)  
VI. [Downloading_the_identified_RNA_Sequences](#Downloading_the_identified_RNA_Sequences)    
VII. [RNA_Sequences_Alignment(Needleman-Wunsch_Algorithm)](#RNA_Sequences_Alignment(Needleman-Wunsch_Algorithm))  
VIII. [Results](#Results)  
IX. [References](#References)  


## I. Introduction:

In this project, I have retrieved and analyzed the RNA sequence data of COVID's prominent strains, Delta and Omicron. 
RNA, a nucleic acid found within every living cell, comprises a single strand composed of diverse arrangements of four nucleotides: uracil, cytosine, adenine, and guanine. 
This RNA serves as the fundamental blueprint for COVID, facilitating cellular entry and self-replication of the virus.

## II. Dataset:

Using the available data from NIH (National Institutes of Health), the federal government agency in the U.S.
Documentation of the data is available here: https://www.ncbi.nlm.nih.gov/datasets/docs/v1/data-packages/sars-cov-2-genome/
The drive link for the dataset can be found here: https://drive.google.com/file/d/1S2ZDjdRkY78kZxBtc9YNUh0mByTHXQ23/view

The Dataset consisted of 847,791 records and  17 attributes.

## III. Data_Preparation:

Firstly, the metadata of Covid RNA Sequences was loaded in the Pandas data frame.
Then, the datatypes of each column in the dataset were checked for appropriateness. "Collection Date" column into the appropriate date datatype. 
And the datatypes of the rest of the columns seemed suitable for analysis.
Further, the "Geolocation" column was split to obtain the continent name and country name in two separate columns.

## IV. Exploratory_Data_Analysis:

1.) At what point in time was the first COVID RNA sequence gathered from each of the continents?

--> 
![image](https://github.com/Priyank0Gandhi/Analyzing-covid-RNA-sequences/assets/96395339/e9c17e60-8bc1-42df-8d89-d8cdcb43112c)

2.) How many sequences were collected on each continent?

--> 
![image](https://github.com/Priyank0Gandhi/Analyzing-covid-RNA-sequences/assets/96395339/28982409-ff2a-4e12-a928-19150e5e4600)            
 ![image](https://github.com/Priyank0Gandhi/Analyzing-covid-RNA-sequences/assets/96395339/9eb2d24a-7f0f-40e0-ac5e-395801b0e1e3)

3.) How long are the shortest and longest sequences? Look at the outliers, if any, to estimate their representation.

--> 
![image](https://github.com/Priyank0Gandhi/Analyzing-covid-RNA-sequences/assets/96395339/ec57a92a-d69a-43df-aebd-8963f8ca8fd5)
    
    Values less than 29639 can be considered outliers.
  Filtered out all the outlier records from the data frame.
    
4.) How many samples were collected by month? Are there any trends?

--> 
![image](https://github.com/Priyank0Gandhi/Analyzing-covid-RNA-sequences/assets/96395339/5c44efe1-4f24-4fde-9a2a-bae0ececf745)

## V. Selecting_RNA_for_analysis:

Now, we'll analyze the sequences themselves rather than the metadata. To do that, we first need to find some sequences we want to analyze.
Let's consider the following sequences:

The reference sequence: the first COVID genome that was fully sequenced.

A base sequence: The first sequence for North America.

Delta sequences (one of the most common COVID variants).

Omicron sequences (another common COVID variant).

I was able to filter out a total of 6 sequences, based on the above conditions.

Pulled out that metadata in a new Pandas Dataframe.

## VI. Downloading_the_identified_RNA_Sequences:

Used the "Biopython" library to retrieve the RNA sequence data through their accession no.

Converted the genomic sequences from FASTA to a more readable format, to be able to work easily during the Analysis.

## VII. RNA_Sequences_Alignment(Needleman-Wunsch_Algorithm):

Aligning the retrieved RNA sequences facilitated analyzing the mutations in RNA sequences over a period of time and in different places of the world.

*Mutation: Alteration of the nucleic acid sequence in the RNA resulting into the formation of a variant*

There are three patterns of Mutations that can occur:/n
i.) Deletion: A nucleotide can be missing from the variant but was present in the original/ reference sequence.
ii.) Substitution: A nucleotide can be replaced by some other sequence of the same size.
iii.) Insertion: A new nucleotide can be added to the variant that is not present in the original/ reference sequence.

We need to determine which parts RNA sequence align and which don't

For this purpose, the "Needleman Wunsch" algorithm was used.

Alignment scores were generated and checked using a 6x6 matrix between all the RNA sequences that were retrieved based on the metadata analysis.

![image](https://github.com/Priyank0Gandhi/Analyzing-covid-RNA-sequences/assets/96395339/54326d60-0a0f-4599-bf10-cb0d2d574a02)

## VIII. Results:

Applied color coding to show the type of mismatch(deletion, substitution, insertion) more clearly.
For this, I used color coding in the HTML format using IPython.

Black- points of sequences that align
Red-   points deleted in the Omicron Variant
Green- points inserted in the Omicron Variant
Blue-  points substituted in the Omicron Variant 
![image](https://github.com/Priyank0Gandhi/Analyzing-covid-RNA-sequences/assets/96395339/d0a41472-2bb5-434e-bd27-939aff8d9da2)

## IX. References:

1.) NIH Datasets: https://www.ncbi.nlm.nih.gov/datasets/docs/v1/reference-docs/data-packages/sars-cov-2-genome/






