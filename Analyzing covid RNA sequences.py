#!/usr/bin/env python
# coding: utf-8

# In[1]:


from platform import python_version
python_version()


# In[2]:


get_ipython().system('pip install ipython')
get_ipython().system('pip install biopython')
get_ipython().system('pip install pandas')


# In[3]:


import pandas as pd
import scipy as sc
import matplotlib.pyplot as plt
import seaborn as sb
import plotly as pt
import numpy as np


# In[4]:


md=pd.read_csv('C:\\Users\\PRIYANK GANDHI\\Downloads\\ncbi_datasets.csv', low_memory=False)
md.head(10)


# In[5]:


md.dtypes


# In[6]:


#Conversion of Date column into appropriate date datatype:
md["Collection Date"]=pd.to_datetime(md["Collection Date"])
md.info()


# In[7]:


#Splitting the GeoLocation to seperate the continent name and country name.
md_geo=pd.DataFrame()
md_geo[["Continent","Country"]]= md["Geo Location"].str.split(";",expand=True)
md_geo["Collection_Date"]=md["Collection Date"]
md_geo


# In[8]:


md_geo.info()


# In[9]:


#EDA1: When was the first COVID RNA sequence collected on each continent?
#grouping by Continents and finding the minimum date
e1=pd.DataFrame()
e1=md_geo.groupby(['Continent']).agg(min_date=('Collection_Date',np.min)).sort_values("min_date")
e1


# In[10]:


plt.plot(e1)


# In[11]:


#EDA2: How many sequences were collected on each continent?

#size()/value_counts() to count the no. of occurances

e2=md_geo['Continent'].value_counts().sort_values()
#md_geo['Continent'].value_counts().sort_values()

#Converting series to a dataframe
e2=pd.DataFrame(e2)

#Continent was considered as an index, so converted the index into a column.
e2=e2.reset_index()

#renaming count column as it does not display any name:
e2.rename(columns={"Continent":"count","index":"Continent"}, inplace=True)
e2


# In[12]:


sb.barplot(x="Continent",y="count", data=e2)


# In[13]:


#EDA3: How long are the shortest and longest sequences? Look at the outliers, if any, to estimate their representation?

md["Nucleotide Length"].describe()

#ANS: Values less than 29639 are considered as outliers[Q1-1.5(IQR)] and more than 30031[Q3+1.5(IQR)]


# In[14]:


#Filtering out all outlier rows
e3=md[(md["Nucleotide Length"]<29639) | (md["Nucleotide Length"]>30031)].sort_values("Nucleotide Length")
e3


# In[15]:


#EDA4: How many samples were collected by month? Are there any trends?

# splitting the date into year, month, and day
e4=md_geo["Collection_Date"].dropna()
e4=pd.DataFrame(e4)
e4["yr","mnth","day"]=e4["Collection_Date"].apply(lambda x: x.timetuple()[:3])
e4.rename(columns={1:"split"}, inplace=True)
#e4.assign(**dict(zip('ghi', e4.col.str)))
e4.head()


# In[16]:


#To generate a continent column-(.str.replace(";.+",'')--This RegEx says replace ; and anything after it with blank'')
md["continent"]=md["Geo Location"].str.replace(";.+",'')
md


# In[17]:


md.groupby("continent").apply(lambda x: x.sort_values('Collection Date').iloc[0])


# In[18]:


#To Identify Sequences to download, generate a dataframe of selected rows.

seq = md[(md['Sequence Type'] == 'RefSeq') | md['Isolate Name'].str.contains('Omicron').fillna(False) | md['Isolate Name'].str.contains('Delta').fillna(False)|(md['Nucleotide Accession'].str.contains('OL467832.1'))]
#Put each condition in parenthesis to ensure correct results
seq


# In[19]:


from Bio import Entrez
Entrez.email="pgwork02@gmail.com"


# In[20]:


#To download sequences from nih database
def download_sequence(id_code):
    handle = Entrez.esearch(db="nucleotide", term=id_code, retmax="1")
    record = Entrez.read(handle)
    handle = Entrez.efetch(db="nucleotide", id=record["IdList"][0], rettype="fasta", retmode="text")
    return handle.read()


# In[41]:


seq_names= seq["Nucleotide Accession"].values.tolist()
seq_names


# In[21]:


seq_data = {}
for sequence in seq["Nucleotide Accession"] :
    seq_data[sequence] = {"fasta": download_sequence(sequence)}


# In[22]:


seq_data


# In[23]:


from Bio import SeqIO
import io

for m,n in seq_data.items():
    f = io.StringIO(n["fasta"])
    seq_data[m]["parsed"] = list(SeqIO.parse(f, "fasta"))[0]


# In[24]:


seq_data['OM108163.1']['parsed']


# In[25]:


#To figure out score to know how different 2 sequences are from each other:

from Bio import Align

aligner= Align.PairwiseAligner()


# In[26]:


#Algorithm name(Needleman-Wunsch): Computer looks at two RNA sequences and see which parts overlap.
aligner.algorithm


# In[27]:


#gives the score of alignment, 1st is the reference sequence
score = aligner.score(seq_data["NC_045512.2"]["parsed"].seq, seq_data["OM061695.1"]["parsed"].seq)
score


# In[33]:


len(seq_data["NC_045512.2"]["parsed"].seq)


# In[34]:


#2nd sequence is aligned 99.71%. to the reference sequence.
29818/29903


# In[45]:


#Check the scores between all the sequences
comparisons = np.zeros((6,6))

for i in range(0,6):
    for j in range(0,i+1):
        score = aligner.score(seq_data[seq_names[i]]["parsed"].seq, seq_data[seq_names[j]]["parsed"].seq)
        comparisons[i,j] = score


# In[46]:


comparisons


# In[50]:


human_names=['reference', 'Delta-3','Delta-1','NorthAmerica','Delta-2','Omicron']


# In[52]:


#Converted numpy array(Matrix) into pandas Dataframe.
comparison_df=pd.DataFrame(comparisons, columns=human_names, index=human_names )


# In[53]:


comparison_df


# In[54]:


#How the sequences align to the reference sequence:
comparison_df.iloc[:,0]/29903


# In[55]:


#Sequence Mismatches(Mutation points)
s1 = seq_data["NC_045512.2"]["parsed"].seq
s2 = seq_data["OM095411.1"]["parsed"].seq
omi_align = aligner.align(s1,  s2)


# In[56]:


omi_align


# In[63]:


omi_aligns=omi_align[0]


# In[64]:


omi_aligns.shape


# In[65]:


omi_aligns.aligned


# In[66]:


#Loop through the 2 sequences and figure out mismatches:
#zip takes 1st element from both the seq and wants them together in one list, then 2nd element from both and puts them in list 2, and so on)
s1_end = None
s2_end = None
for alignments in zip(omi_aligns.aligned[0], omi_aligns.aligned[1]):
    
    if s1_end and s2_end:
        s1_mismatch = s1[s1_end:alignments[0][0]]
        s2_mismatch = s2[s2_end:alignments[1][0]]
        print("1: {}".format(s1_mismatch))
        print("2: {}".format(s2_mismatch))
    
    s1_end = alignments[0][1]
    s2_end = alignments[1][1]


# In[67]:


from IPython.display import HTML


# In[68]:



def color_print(s, color='black'):
    return "<span style='color:{}'>{}</span>".format(color, s)


# In[69]:


s1_end = None
s2_end = None
display_seq = []
for alignments in zip(omi_aligns.aligned[0], omi_aligns.aligned[1]):
    
    if s1_end and s2_end:
        s1_mismatch = s1[s1_end:alignments[0][0]]
        s2_mismatch = s2[s2_end:alignments[1][0]]
        if len(s2_mismatch)==0:
            display_seq.append(color_print(s1[s1_end:alignments[0][0]], "red"))
        elif len(s1_mismatch)==0:
            display_seq.append(color_print(s2[s2_end:alignments[1][0]], "green"))
        else:
            display_seq.append(color_print(s2[s2_end:alignments[1][0]], "blue"))
    display_seq.append(s1[alignments[0][0]:alignments[0][1]])
    s1_end = alignments[0][1]
    s2_end = alignments[1][1]


# In[70]:


display_seq = [str(i) for i in display_seq]


# In[71]:


display(HTML('<br>'.join(display_seq)))


# In[ ]:




