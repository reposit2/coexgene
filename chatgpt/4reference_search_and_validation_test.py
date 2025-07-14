# Notebook for running literature search and validation

import os
import pandas as pd
import json 
from Bio import Entrez
import openai
from utils.reference_checker import get_references_for_paragraphs, iter_dataframe
import pickle
#%load_ext autoreload
#%autoreload 2

runVersion = 'initial'#'initial'
#runVersion = 'additional'

dataType = "Omics"

if dataType == "Omics":
    LLM_analysisFilePath = 'data/omics_revamped_LLM_DF_test.tsv'
    toSaveFilePath = 'data/omics_revamped_LLM_ref_DF_test' # remove the .tsv, because output will be saved as a dataframe(.tsv) and a dictionary (.json) in the function
    jsonFilePath = 'jsonFiles/reference_checking_revision_test.json'
    nameCol = 'GeneSetName'
    LLM_analysisCol = 'LLM Analysis'

elif dataType == "GO_sets":
    LLM_analysisFilePath = 'data/GO_term_analysis/simrank_LLM_processed_selected_1000_go_terms.tsv'
    toSaveFilePath = 'data/GO_term_analysis/simrank_LLM_processed_selected_1000_go_terms_refs'
    jsonFilePath = 'reference_checking_task1.json'
    examplesTORun = ["GO:0019433"] 
    nameCol = 'GO'
    
else:
    print("Not implemented for dataType")

if runVersion == 'additional':
    LLM_analysisFilePath = toSaveFilePath + '.tsv'

LLM_analysisFilePath

# OUT [9]: 'data/omics_revamped_LLM_DF.tsv'

with open(jsonFilePath) as json_file:
    config = json.load(json_file)

openai.api_key = os.environ["OPENAI_API_KEY"]
email = config['EMAIL']

## Read in the LLM analysis file
df = pd.read_csv(LLM_analysisFilePath, sep='\t', keep_default_na=False, na_values=['NaN'])
df = df.replace({'None':None})
df.head()

# OUT [11]:
#Source	GeneSetID	GeneSetName	GeneList	n_Genes	LLM Name	LLM Analysis	Score
#0	NeST	Cluster1-10	Cluster1-10	CTRL HSD17B14 KIAA0232 PAQR8 PLA2G1B RNF145 SG...	12	Lipid Metabolism and Membrane Dynamics	1. CTRL, or carboxyl-terminal esterase/lipase,...	0.85
#1	NeST	Cluster1-11	Cluster1-11	LMF1 MFHAS1 MR1 PLA2G1B RASL11A RNF145 SLC2A6 ...	12	Lipid Metabolism and Membrane Trafficking	1. LMF1 (Lipase Maturation Factor 1) is crucia...	0.85
#2	NeST	Cluster1-12	Cluster1-12	AMY2B CNPY2 EGFL7 LDLR LPL LRP8 LRPAP1 MYLIP P...	12	Lipid Metabolism and Receptor-Mediated Endocyt...	1. AMY2B (Amylase, alpha 2B) is an enzyme that...	0.88
#3	NeST	Cluster1-13	Cluster1-13	ACO1 CDKAL1 CDX2 CIAO1 CITED2 FAM96A FAM96B FB...	16	Iron Metabolism and Regulation in Pancreatic F...	1. ACO1, also known as aconitase 1, is a key e...	0.88
#4	NeST	Cluster1-14	Cluster1-14	DTWD1 DTWD2 POLD2 POLD4 POLH POLI RAD18 REV1 R...	15	DNA Damage Tolerance and Repair Pathway	1. DTWD1 and DTWD2 are proteins containing DTW...	0.92

## UPDATES 02/13/2023
#iterate through the df and save df and dict

res_df, res_dict = iter_dataframe(df, email, config, n=3, papers_query=20, verbose=False, return_paragraph_ref_data=False, id_col=nameCol, paragraph_col=LLM_analysisCol, runVersion=runVersion, save_path = toSaveFilePath)

## merge the referenced analysis to the reduced_LLM_genes_APV_only_DF
#reduced_LLM_genes_APV_only_DF = pd.read_csv('data/omics_revamped_LLM_gprofiler_new_gene_name_DF_APV_only_test.tsv', sep="\t")
#referenced_analysis_file = 'data/omics_revamped_LLM_ref_DF_test.tsv'
#referenced_analysis_df = pd.read_csv(referenced_analysis_file, sep="\t")[['Source',
#       'GeneSetID', 'GeneList', 'n_Genes','LLM Name', 'referenced_analysis']]
## referenced_analysis_df.shape
## referenced_analysis_df.columns
## make sure the values geneset ID and LLM name are the same with the reduced_LLM_genes_APV_only_DF
#reduced_LLM_genes_APV_only_DF.loc[(reduced_LLM_genes_APV_only_DF['GeneSetID'] == referenced_analysis_df['GeneSetID']) & (reduced_LLM_genes_APV_only_DF['LLM Name'] == referenced_analysis_df['LLM Name']), ['referenced_analysis']] = referenced_analysis_df['referenced_analysis']
#reduced_LLM_genes_APV_only_DF[['GeneSetID', 'GeneList', 'n_Genes','LLM Name', 'LLM Analysis', 'referenced_analysis']]

# OUT [18]:
#GeneSetID	GeneList	n_Genes	LLM Name	LLM Analysis	referenced_analysis
#0	BRD-A00546892_-666_MCF7_6.0_h_10.0_um	CITED2 COL5A1 CRABP2 KCTD12 MDFIC MMP2 NRP1 OR...	24	Cellular Matrix Remodeling and Tissue Development	1. CITED2, TWIST1, and LMO2 are transcriptiona...	1. CITED2, TWIST1, and LMO2 are transcriptiona...
#1	BRD-A00993607_ALPRENOLOL_MCF7_6.0_h_10.0_um	1060P11.3 ADM AHR AMIGO2 ARL4C ATP10D CAV2 CD4...	47	Cellular Adhesion and Extracellular Matrix Int...	1. Several proteins in this system, such as CD...	1. Several proteins in this system, such as CD...
#2	BRD-A00993607_Alprenolol hydrochloride_MCF7_6....	ABAT ASS1 CHI3L1 CHST2 CLDN3 EIF5B FRZB GAL HE...	59	System of unrelated proteins	The provided list of interacting proteins enco...	The provided list of interacting proteins enco...
#3	BRD-A01320529_Salmeterol_MCF7_6.0_h_10.0_um	AMIGO2 AREG GAS6 GPR37 IFT57 PELI1 SQLE AKAP12...	9	System of unrelated proteins	1. AMIGO2 (Amphoterin Induced Gene and ORF) is...	1. AMIGO2 (Amphoterin Induced Gene and ORF) is...
#4	BRD-A01346607_FLUMETHASONE_MCF7_6.0_h_10.0_um	1060P11.3 CPE EFNB2 HIST1H2AC IL1R2 INHBB LYPD...	15	System of unrelated proteins	1. The proteins listed do not appear to conver...	1. The proteins listed do not appear to conver...
#...	...	...	...	...	...	...
#295	hMPV_72Hour	RNMT SAMD9L CEBPB CEBPG PSMD12 RSRC2 XBP1 GADD...	22	Cellular Stress Response and Protein Homeostasis	1. RNMT, the RNA (guanine-7-) methyltransferas...	1. RNMT, the RNA (guanine-7-) methyltransferas...
#296	icSARA deltaORF6_48Hour	FGF19 ADM2 LRIT1 UCHL1 SLC19A1 TSSC1 INPP5J HY...	49	System of unrelated proteins	The provided list of interacting proteins enco...	The provided list of interacting proteins enco...
#297	icSARS CoV_12Hour...321	TSPY3 LGALS14 PPBPP2 LARS MRPS10 PAX9 RBM22 LO...	16	System of unrelated proteins	1. TSPY3, testis-specific protein Y-linked 3, ...	1. TSPY3, testis-specific protein Y-linked 3, ...
#298	icSARS CoV_3Hour	ANK2 GAB3 TGFB3 CECR3 MUC19 LOC100507053 ZNF36...	99	System of unrelated proteins	The provided list of interacting proteins enco...	The provided list of interacting proteins enco...
#299	icSARS CoV_72Hour...15	KRT75 ZDBF2 HIF3A TFAP2D F12 CASP14 FAM151B NG...	72	System of unrelated proteins	The provided list of interacting proteins enco...	The provided list of interacting proteins enco...

# 300 rows × 6 columns

# use your own iterate dataframe loop

#LLM_analysisFilePath = './data/omics_revamped_LLM_DF.tsv'
#toSaveFilePath  = './data/test_omics_revamped_LLM_DF_refs'
#jsonFilePath = 'jsonFiles/reference_checking_revision_test.json'
##with open(jsonFilePath) as json_file:
#    config = json.load(json_file)
#
#openai.api_key = os.environ["OPENAI_API_KEY"]
#email = config['EMAIL']
#nameCol = 'GeneSetName'
#LLM_analysisCol = 'LLM Analysis'
#runVersion = 'additional'

#runOnlyExamples = True
#if runOnlyExamples:
#    examplesTORun = ['BRD-A00993607 Alprenolol hydrochloride MCF7 6.0 h 10.0 um',
#       'BRD-A13964793 -666 MCF7 6.0 h 10.0 um',
#       'BRD-A19633847 PERHEXILINE MALEATE MCF7 6.0 h 10.0 um',
#       'BRD-A31204924 -666 MCF7 6.0 h 10.0 um', 'Cluster2-126',
#       'Cluster2-140', 'Cluster2-147', 'Cluster2-169', 'Cluster2-183',
#       'Cluster2-191', 'Cluster2-200']

## Read in the LLM analysis file
df = pd.read_csv(LLM_analysisFilePath, sep='\t', keep_default_na=False, na_values=['NaN'])
df = df.replace({'None':None})
df.set_index(nameCol, inplace=True)
# # rename 'References' to 'referenced_analysis' 
# df = df.rename(columns={'References':'referenced_analysis'})

# Out [11]:
#	Unnamed: 0.1	...1	...2	Unnamed: 0	Source	GeneSetID	GeneList	n_Genes	LLM Name	LLM Analysis	...	Term	GO term	GO ID	GO_term_genes	LLM_name_GO_term_sim	enrichr_JI	LLM_success_TF	enrichr_success_TF	referenced_analysis	enrichr_success_TF_0.1
#GeneSetName																					
#BRD-A31204924 -666 MCF7 6.0 h 10.0 um	85	85	9762	1	L1000	BRD-A31204924_-666_MCF7_6.0_h_10.0_um	1060P11.3 ANPEP DPP4 FAM129A HNMT MKNK1 MUC5B ...	15	System of unrelated proteins	1. 1060P11.3 is a hypothetical protein with li...	...	Histamine Metabolic Process (GO:0001692)	Histamine Metabolic Process	GO:0001692	SLC22A3 PRG3 TRH HNMT SLC29A4	0.264298	0.052632	False	False		False
#1 rows × 26 columns

#if runVersion == 'initial':
#    df['referenced_analysis'] = None

#####USE get_references_for_paragraph####
#saved_dict = {}
#i = 0
#for set_id, row in df.iterrows():
##for i in range(startRow, df.shape[0]):
##    row = df.iloc[i]
#    # if runOnlyExamples: # Only run examples
#    #    if df.iloc[i][nameCol] not in examplesTORun: 
#    #        continue
#    if runVersion == "initial":
#        if df.loc[set_id, 'referenced_analysis'] is not None:
#            continue
#    if runVersion == "additional":
#        with open(toSaveFilePath + '.json') as f:
#            saved_dict = json.load(f)
#        if not (df.loc[set_id,'referenced_analysis'] == ''):
#            continue # skip this row because already done
#            
#    print('=========================================')
#    print('=========================================')
#    print('=========================================')
#
#    print(['dataframe row', set_id])
#    # check out the llm analysis 
#    example_analysis = df.loc[set_id, LLM_analysisCol]
#    paragraphs = list(filter(lambda p: len(p.split()) > 5, example_analysis.split("\n")))
#    
#    try:
#        references, paragraph_dict = get_references_for_paragraphs(paragraphs, email = email, config =config, n=3, verbose=True, papers_query=20, return_paragraph_ref_data=True)
#        
#    except Exception as e:
#        print('Cannot get references for row', set_id, e)
#        references = ''
#        saved_dict[set_id] = None
#        
#    references.replace( '\n', '')
#    
#    df.loc[set_id, 'referenced_analysis'] = references
#    saved_dict[set_id] = paragraph_dict
#    i += 1
#    if i%5==0:
#        df.to_csv(toSaveFilePath+ '.tsv', sep = '\t')
#        with open(toSaveFilePath + '.json', 'w') as f:
#            json.dump(saved_dict, f)
#        
## if not runOnlyExamples: 
#df.to_csv(toSaveFilePath+'.tsv', sep = '\t')
#with open(toSaveFilePath + '.json', 'w') as f:
#    json.dump(saved_dict, f)
## check there is no None
#print(len(df[df['referenced_analysis'] == '']))
#
#df.to_csv(toSaveFilePath+'.tsv', sep = '\t')

