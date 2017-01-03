#Importing the necessary libraries
import argparse
import os
import pandas as pd
import datetime


#### Constants for the pipeline
threshHLA = 0.75

#### Recommendation databases
drugVatiants = pd.read_csv('Data/total_var_ann.csv')
HLAVariants = pd.read_csv('Data/HLA_var_ann.csv')

#Implementing arguments importation from the command line type
#python main.py -f filname to import filename into the programm.
parser = argparse.ArgumentParser(description='Import the file name into the programm')
parser.add_argument("-f", "--file", dest="filename",
                  help="Option to add filename", metavar="FILE")

args = parser.parse_args()

#Define the patient file
patientFile = args.filename
patientName = patientFile.split('.')[0]

#Make a folder for the results of the patient
#We store the reults in separate folders each day there is an analysis
today = str(datetime.datetime.today().date())
path_today = 'Results/'+today
if not os.path.exists(path_today):
    os.mkdir(path_today)

#We also create a directory for the patient name
pathPatient = path_today+'/'+patientName
if not os.path.exists(pathPatient):
    os.mkdir(pathPatient)

#Converting the file of the patient to plink
os.system('plink --vcf '+patientFile+' --recode bimbam --out '+patientName+' --make-bed --double-id')

#Extracting the chromosome 6
patientChrom6 = patientName+'_chr6'
os.system('plink --bfile '+patientName+' --chr 6 --out '+ patientChrom6 +' --make-bed')

#COmputing HLA imputation with HIBAG method
os.system('Rscript --vanilla HLA_from_patient.R '+patientChrom6+' '+pathPatient)

#move HLA results to the right folder.
hlaPatient = pd.read_csv(patientChrom6+'_results.csv')
os.system('mv '+patientChrom6+'_results.csv '+pathPatient)

#Now we import tne results file so that we modify it to be more readable
hlaPatient.columns = hlaPatient.loc[0]
hlaPatient.drop(0,inplace=True)
HLAFinal = hlaPatient.copy()
hlaPatient.drop(1,axis=1,inplace=True)
hlaPatient.to_csv(pathPatient+'/ResultsHLA.csv')
#Now we filter the locuses to be sure that all locuses are rightly imupted

for locus in hlaPatient.columns:
    if float(HLAFinal.loc[3,locus])<threshHLA:
        HLAFinal.drop(locus,axis=1,inplace=True)

#Modifying the datafame so that we can
for locus in HLAFinal.columns[1:]:
    allele1 = HLAFinal.loc[1, locus].split(':')
    allele2 = HLAFinal.loc[2, locus].split(':')

    HLAFinal.loc[1, locus] = 'HLA_' + locus + '_' + allele1[0] + allele1[1]
    HLAFinal.loc[2, locus] = 'HLA_' + locus + '_' + allele2[0] + allele2[1]

#Matching with the HLA database
arr_HLA = HLAVariants['Variant']
HLAPAtient = []
for index in HLAVariants.index:
    locus = HLAVariants.loc[index,'Variant']
    for col in HLAFinal.columns[1:]:
        if HLAFinal[col].isin([locus]).any()==True:
            HLAPAtient.append(HLAVariants[HLAVariants.index==index])

PatientHLA = pd.concat(HLAPAtient)
PatientHLA = PatientHLA[['Variant','Chemical','Note']]

#Saving the result in a csv file
PatientHLA.to_csv(pathPatient+'/Results_HLA_Recommendations.csv')

os.system('mv '+patientFile+' '+pathPatient)

#Clean the directory
for file in os.listdir():
    if file.startswith(patientName):
        os.remove(file)

