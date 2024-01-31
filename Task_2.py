#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
df = pd.read_csv('skin_toxicity.csv')
df.head(10)


# **Molecular descriptor calculation using RDKit**

# In[2]:


from rdkit import Chem
from rdkit.Chem import Descriptors


# In[3]:


# Define a function to calculate descriptors from SMILES
def calculate_all_descriptors(smiles):
    # Convert SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    # If the molecule is invalid, return None for all descriptors
    if mol is None:
        return {desc_name: None for desc_name in Descriptors._descList}
    # Otherwise, calculate descriptors and store them in a dictionary
    descriptors = {}
    for desc_name, desc_func in Descriptors._descList:
        try:
            descriptors[desc_name] = desc_func(mol)
        except:
            descriptors[desc_name] = None
    # Return the dictionary of descriptors
    return descriptors


# In[4]:


# Apply the function to the SMILES column
df_descriptors = df['SMILES'].apply(calculate_all_descriptors)


# In[5]:


# Convert the list of dictionaries to a data frame
df_descriptors = pd.DataFrame(df_descriptors.tolist())


# In[6]:


# Merge with the original data frame
df_combined = pd.concat([df, df_descriptors], axis=1)


# In[7]:


# Print the first 10 rows
print(df_combined.head(10))


# In[8]:


# Save the data frame as excel file
df_combined.to_excel('skin_toxicity_with_descriptors.xlsx', index=False)


# **Using the SiRMS**

# In[9]:


# Import RDKit and pandas
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess


# In[10]:


# Use the canonical_smiles column instead of the SMILES column
df['SMILES'].to_csv('smiles_input.txt', index=False, header=False)


# In[12]:


# Convert SMILES to SDF format
writer = Chem.SDWriter('molecules.sdf')
for smile in df['SMILES']:
    mol = Chem.MolFromSmiles(smile)
    if mol is not None:
        AllChem.Compute2DCoords(mol)
        writer.write(mol)
writer.close()


# In[13]:


# Define the SiRMS command
command = 'sirms -i molecules.sdf -o output02.txt --min_atoms 2 --max_atoms 6 -x'

# Run the SiRMS command
process = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


# In[14]:


# Read the output file
output_df = pd.read_csv('output02.txt', sep='\t') 


# In[ ]:


# Print the first 100 rows
output_df.head(100)


# In[ ]:


# Print the information of the output data frame
output_df.info()


# In[15]:


# Save the SMILES column as csv file for OCHEM
df['SMILES'].to_csv('smiles_for_ochem.csv', index=False, header=False)


# **I caculate rest of chosen 2D descriptors using Ochem**

# In[16]:


# Load the OChem calculated descriptors
df_ochem = pd.read_csv('ochem-descriptors.csv')
df_ochem.head(10)


# In[17]:


# merge descriptor calculation using RDKit with rest of chosen 2D descriptors using Ochem
df_combined_merged = pd.merge(df_combined, df_ochem, on='SMILES', how='left')

df_combined_merged.head(10)


# In[18]:


#to merge descriptor calculation Using the SiRMS first we need common column 
smiles_df = pd.read_csv('smiles_for_ochem.csv', header=None, names=['SMILES'])
output_df['SMILES'] = smiles_df['SMILES']
print(output_df.head())


# In[19]:


final_df = pd.merge(df_combined_merged, output_df, on='SMILES', how='left')

final_df.head(10)


# In[29]:


# Assuming final_df is your DataFrame
final_df = final_df.loc[:, (final_df != 0).any(axis=0)]

# Assuming final_df is your DataFrame with duplicate columns
final_df = final_df.loc[:, ~final_df.columns.duplicated()]
# This will keep only the first occurrence of each column name
# Assuming final_df is your DataFrame with missing or NA values
final_df = final_df.dropna(axis=1, how='any')
# This will drop any column that has at least one missing or NA value



final_df.head(10)


# In[30]:


final_df.to_excel('skin_toxicity_final_descriptors.xlsx', index=False)


# In[ ]:




