from Bio import Phylo
import ete3
import sys, os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from io import StringIO
from matplotlib.colors import ListedColormap




#########################################################
#                   Main Program                        #
#########################################################



# Common variables
specie = []

#miRNA treatment variables
finalData_miRNA = []
start_pos_miRNA = []
miRNA_start_end = []
#DNA sequence treatment variables
finalData_DNA = []
DNA_start_end = []
shift = []
#Newick file reading and tree contruction
modified_names = [] #filtered specie names. The format of the species names variations in the newick files : "familyname_speciename_XXXX" or "familyname_speciename_XXXX_YYYY" or "familyname_speciename_XXXX_YYY_ZZZ". What whe want : familyname_speciename


filename = sys.argv[1]
newick_format = sys.argv[2]

##########################################################################################################################################
#                                                       Dataframes Creation                                                              #
##########################################################################################################################################


with open(filename,"r") as file :
    lines = file.readlines()
    print("###############################################")
    print("I am reading the file")
    for line in lines :
        ################################################################################
        #extract species name : needed for both dataframes
        specie_found = line.split("\t")[0]
        specie.append(specie_found)

        ################################################################################
        #tmp contains all the rest of the line excluding the specie name (aka matches, starting_pos and end_pos)
        tmp = line.split("\t")[1]
        Rawdata = tmp.split()[0]
        #we extract the matches of the miRNA in the file
        finalData_found_miRNA = Rawdata.split("&")[0]
        finalData_miRNA.append(finalData_found_miRNA)

        #Now we do the same but for the DNA sequence matches, both are just separated by an "&"
        finalData_found_DNA = Rawdata.split("&")[-1]


        ################################################################################
        #extract the start_pos of the miRNA
        miRNA_start_extract = tmp.split(":")[0]
        miRNA_start_end_pos = miRNA_start_extract.split()[-1]
        #print(miRNA_start_end_pos)
        miRNA_start_end.append(miRNA_start_end_pos)

        miRNA_start_pos = miRNA_start_end_pos.split(',')[0]
        start_pos_miRNA.append(miRNA_start_pos)


        ################################################################################
        #extract the start_pos on the DNA_sequence, used to get the number between start_pos and first nucleotide of miRNA : to fill the stored sequences for the heatmap so that all sequences are aligned and compared from the same nucleotide
        DNA_start_extract = tmp.split(":")[-1]
        DNA_start_end_pos = DNA_start_extract.split()[0]
        DNA_start_end.append(DNA_start_end_pos)
        DNA_end_pos = DNA_start_end_pos.split(',')[-1] #extract the end position, needed to get the number between the end_pos for the ORTHOLOGOUS and the lenght of the sequence so we can get the shift value required for the alignment of the heatmap
        ################################################################################

        #Calculate the shift we have to make, this operation is for visual purpose.It concerns the orthologous sequence, where we want the matching nucleotide with the miRNA seed to be aligned on the heatmap 
        starting_pos_DNA= DNA_start_extract.split(",")[0]

        # to get the shift value required for the coherence of the heatmap : if the shift is positive we make a left shift by deleteing N-values on the right of the sequence, else if the shift is negative we make a right shift by adding N-values
        if (miRNA_start_pos =="null") :
            shift_for_heatmap = "0" 
        else :
            if (miRNA_start_pos =="1") :
                    if int(DNA_end_pos) >= 40 :
                        shift_for_heatmap = int(DNA_end_pos) - 40 - int(miRNA_start_pos)
                    else :    
                        shift_for_heatmap = int(DNA_end_pos) - 40 - int(miRNA_start_pos)
            else :
                if int(DNA_end_pos) >= 40 :
                    shift_for_heatmap = int(DNA_end_pos) - 40 + int(miRNA_start_pos) -2
                else :
                    shift_for_heatmap = int(DNA_end_pos) - 40 + int(miRNA_start_pos) -2 


        #print("l'espèce : ",specie_found,"\nla position start_end de l'ADN : ",DNA_start_end_pos,"\nla ed_position : ",DNA_end_pos,"\n le shift :",shift_for_heatmap,"\n\n")
        #print("le calcul : end_pos - anchor pos(40) : ",DNA_end_pos," -40 ; + miRNA_start_pos -2 : ",miRNA_start_pos," -2"," = shift : ",shift_for_heatmap,"\n\n")

        shift.append(shift_for_heatmap)

        ################################################################################
        # Now we normalize the data sont that on the heatmap every sequence starts at the same point as the others (some miRNA are aligned on nucleotide 2 of the DNA sequence, others on 3,4,7,16 ... etc) 
        if starting_pos_DNA == "null\n" :
            shiftWide = "0"
        else : 
            shiftWide = int(starting_pos_DNA)-1   
        #We add as much characters to fill the gaps, so that if a miRNA is matching with nucleotide 1 : no char is added. And if it matches the N-nucleotide we add N-1 characters
        finalData_found_DNA = ((int(shiftWide) - int(shift_for_heatmap) +1) * "#") + finalData_found_DNA 
        finalData_DNA.append(finalData_found_DNA)

        ################################################################################

print("\ndata extracted with success")        
print("\nShift values calculated with success")
print("###############################################")



##############################################################  miRNA DataFrame Creation   ##########################################################################################

df = pd.DataFrame(columns =["species","matches_miRNA"])        
df["species"]= specie
df["matches_miRNA"]= finalData_miRNA
df["position"] = miRNA_start_end
df["start_pos_miRNA"] = start_pos_miRNA

pd.set_option('display.max_rows', df.shape[0]+1)
#print(df)


###### Here we create a sub dataframe that contains the 4 species that have duplicates sequences. Since we are using a dictionnary later on, the duplicates are 
# deleted since the keys are unique. A simple strategy is to merge this dataframe with the final one (used to make the heatmap after all the filtering steps) 
# and give every duplicate_specie the same index as the original specie. ######

duplicate_species = df[df['species'].str.endswith('_2')]



######### this part is to eliminate species with "null" values for an alternative heatmap : optional #########

#df = df[df["matches_miRNA"].str.contains("null")==False]
#df = df.reset_index()


pd.set_option('display.max_rows', df.shape[0]+1)
#print(df)
##############################################################  DNA sequence DataFrame Creation   ##########################################################################################




df_seq = pd.DataFrame(columns =["species","matches_DNA","position"])        
df_seq["species"]= specie
df_seq["matches_DNA"]= finalData_DNA
df_seq["position"]= DNA_start_end
df_seq["shift"]= shift
duplicate_seq_species = df_seq[df_seq['species'].str.endswith('_2')]


#print(df_seq)

print("Both dataframes were created")
##########################################################################################################################################
#                                           Newick Data input and data reorganisation                                                    #
##########################################################################################################################################



with open(newick_format, 'r') as tree_file:
    newick_data = tree_file.read().strip()
    #print(newick_data)

# Parse the Newick data to create the tree object
tree = Phylo.read(StringIO(newick_data), 'newick')



species_tree = ete3.Tree(newick_data, format=1, quoted_node_names=True)
print(species_tree)
species_labels = []
species_order = []
cpt = 1

# Traverse the species tree and extract labels and order
for leaf in species_tree.iter_leaves():
    species_labels.append(leaf.name)
    #print(species_labels)
    
    species_order.append(cpt)
    cpt = cpt+1
    #print(len(species_order))


#########   Here we have two choices to filter the species names depending on the newick file we use : 
# - If the newick file names are written as XXXX_familyName_SpecieName  or FamilyName_SpecieName_XXXX/FamilyName_SpecieName_XXXX_YYY ..Etc we use option 1
# - If the newick file is written as : FamilyName.SpecieName we use option 2
# - If newick file is written as : FamilyName.SpecieName_XXXX_...etc we use option 2 to change every dot to an underscrore AND also option 1 to filter the names 
#########


# Option 1 :
species_labels = [element.split('_', 1)[1] for element in species_labels]

# Option 2 :
#species_labels = [element.replace(".","_") for element in species_labels]
#print(species_labels)



for name in species_labels:
    new_name = ""
    first_underscore_found = False

    for char in name:
        if char == "_":
            if first_underscore_found:
                break
            else:
                first_underscore_found = True
        new_name += char

    modified_names.append(new_name)

#print(modified_names)



common_species = []
specie_found = 0
for x in specie :
    for y in modified_names :
        if x == y :
            specie_found = 1
            common_species.append(x)
    
if specie_found == 0 :
    print("No species found")
else :
    print("found species !\n")
    print(common_species)    



print("\n\n",len(modified_names),"Species found in the newick file !!\n")

print((int(len(common_species))-1)," common species between the newick file and our 95 species !!!")



###############################  Reorganize the order of the miRNA dataframe so that species keep their phylogenetical relations  ###############################


species_order_dict = dict(zip(species_order, modified_names))

#We create a dataframe of the species we read from the tree, and add their index as a column, this column will be usefull to reorganise our first dataframe in a way to keep the phylogenetical relations between the species

df_ranking = pd.DataFrame.from_dict(species_order_dict, orient='index', columns=['ordered_species'])
df_ranking["index"] = df_ranking.index


#Merge the two DataFrames based on the common species column, species is from df(99 species), ordered_species from df_ranking(522 species), but we only keep the common species
merged_df = df.merge(df_ranking, left_on="species", right_on="ordered_species", how="inner")

merged_df = pd.concat([merged_df, duplicate_species], axis=0)


# This line is purely esthetic, so that we put out specie of interest in the top of our heatmap, and the duplicates are near the original species
merged_df.loc[merged_df["species"] == "Oryzias_latipes", "index"] = 1 # to force it to be the first of the list
merged_df.loc[merged_df["species"] == "Lepisosteus_oculatus", "index"] = 9999 # to force it to be the last of the list
merged_df.loc[merged_df["species"] == "Simochromis_diagramma_2", "index"] = 13
merged_df.loc[merged_df["species"] == "Myxocyprinus_asiaticus_2", "index"] = 88
merged_df.loc[merged_df["species"] == "Oncorhynchus_nerka_2", "index"] = 72
merged_df.loc[merged_df["species"] == "Dunckerocampus_dactyliophorus_2", "index"] = 68

#concatenating the two dataframes changed my "index" column from int to float, here I change the value type to int again since it is needed later on
merged_df["index"] = merged_df["index"].astype(int)



#Sort the merged DataFrame based on the order in the second DataFrame
merged_df.sort_values(by="index", inplace=True)
#print(merged_df)


df = merged_df[["species", "matches_miRNA", "position", "start_pos_miRNA"]]



###### Some species are present 2 times (actually 4 species) because they have two sequences that aligns well with the miRNA. Here we add these species to the 
# final dataframe since they were deleted when creating the dictionnary ######
df = df.drop_duplicates(subset='species')



###############################  Reorganize the order of the DNA dataframe so that species keep their phylogenetical relations  ###################################

merged_df_seq = df_seq.merge(df_ranking, left_on="species", right_on="ordered_species", how="inner")

merged_df_seq = pd.concat([merged_df_seq, duplicate_seq_species], axis=0)
print(merged_df_seq)

# This line is purely esthetic, so that we put out specie of interest in the top of our heatmap, and the duplicates are near the original species
merged_df_seq.loc[merged_df_seq["species"] == "Oryzias_latipes", "index"] = 1
merged_df_seq.loc[merged_df_seq["species"] == "Lepisosteus_oculatus", "index"] = 9999
merged_df_seq.loc[merged_df_seq["species"] == "Simochromis_diagramma_2", "index"] = 13
merged_df_seq.loc[merged_df_seq["species"] == "Myxocyprinus_asiaticus_2", "index"] = 88
merged_df_seq.loc[merged_df_seq["species"] == "Oncorhynchus_nerka_2", "index"] = 72
merged_df_seq.loc[merged_df_seq["species"] == "Dunckerocampus_dactyliophorus_2", "index"] = 68

merged_df_seq["index"] = merged_df_seq["index"].astype(int)

merged_df_seq.sort_values(by="index", inplace=True)

df_seq = merged_df_seq[["species", "matches_DNA", "position", "shift","index"]]

df_seq = df_seq.drop_duplicates(subset='species')


# Reset the index of the dataframes
df.reset_index(drop=True, inplace=True)
df_seq.reset_index(drop=True, inplace=True)

print(df_seq)

pd.set_option('display.max_rows', df.shape[0]+1)
#print(df_seq)

print("\nBoth dataframes were reorganized based on the newick file order")
print("###############################################")



##########################################################################################################################################
#                Optional part : to create a tree of our 95 based on the newick file we have (Visualisation purpose)                     #
##########################################################################################################################################





##########################################################################################################################################
#                                                         HEATMAP Creation                                                               #
##########################################################################################################################################

#The script creates two heatmaps : one for the miRNA matches obtained from RNAduplex, we can see for example to which degree the seed sequence is conserved
#between the species.
#The second heatmap in for the DNA sequences that matched with these miRNAs, allowing us to see to which degree the sequences are conserved, and for example 
#if the sequences are more conserved if they match a seed,and also in which clades these sequences are conserved.

##############################################################   miRNA sequence treatment part   ##########################################################################################


# Find the maximum length of sequences
max_length = max(len(seq) for seq in df['matches_miRNA'])

# fill sequences with a "." character (considered as blanc on the heatmap) so that all our sequences are aligned on the index of the nucleotide they match with.
filled_sequences = [seq.ljust(max_length, '.') for seq in df['matches_miRNA']]
#print(filled_sequences)



# Create a color map
mycolors = ListedColormap(['white','black'])
# Create a binary matrix
binary_matrix = np.array([[1 if char == '(' else 0 for char in seq] for seq in filled_sequences])

# Set the figure size
plt.figure(figsize=(400, 300))


# Plot the heatmap
ax = sns.heatmap(binary_matrix, cmap=mycolors, cbar=False)


ax.set_xticklabels(np.arange(1, max_length + 1))
ax.set_xticks(np.arange(0.5, max_length))
# Set the species labels on the left side
ax.set_yticks(np.arange(len(df)))
ax.set_yticklabels(df['species'], rotation=0)
plt.savefig('heatmap_miRNA_matches.svg', format='svg')


#######################################################   DNA sequence treatment part   ##########################################################################################


# Find the maximum length of sequences
max_length2 = max(len(seq) for seq in df_seq['matches_DNA'])
#print(max_length2)

# fill sequences with a "." character (considered as blanc on the heatmap) so that all our sequences are aligned on the index of the nucleotide they match with.
filled_sequences2 = [seq.ljust(max_length2, '.') for seq in df_seq['matches_DNA']]
#print(filled_sequences)




# Create a color map
mycolors = ListedColormap(['white','black'])

# Create a binary matrix
binary_matrix2 = np.array([[1 if char == ')' else 0 for char in seq] for seq in filled_sequences2])

# Set the figure size
plt.figure(figsize=(100, 150))

# Plot the heatmap
ax = sns.heatmap(binary_matrix2, cmap=mycolors, cbar=False,)


ax.set_xticklabels(np.arange(1, max_length2 + 1))
ax.set_xticks(np.arange(0.5, max_length2))
#Set the species labels on the left side
ax.set_yticks(np.arange(len(df_seq)))
ax.set_yticklabels(df_seq['species'], rotation=0)

# Better display of the heatmap
sns.set(font_scale=0.5)

# Show both of the plots 
#plt.show()
plt.savefig('heatmap_DNA_matches.svg', format='svg')
print("\nheatmap created with success")

