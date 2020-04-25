import os, os.path
import pandas as pd
from os import listdir
from os.path import isfile, join
cwd=str(os.getcwd())
init_dir=str(os.getcwd())
os.makedirs(str(cwd)+str("/HEPTAD_SITE_WITH_SCORE/"), exist_ok=True)
os.makedirs(str(cwd)+str("/HEPTAD_SITE_WITH_SCORES_ALL_WINDOWS_GREATER_THAN_10.5/"), exist_ok=True)
os.makedirs(str(cwd)+str("/HEPTAD_SITE_WITH_SCORES_EACH ANTIGENS_BEST_HEPTAD_SCORING_WINDOW/"), exist_ok=True)
os.makedirs(str(cwd)+str("/PERCENT_IDENTITY_PANDAS/"), exist_ok=True)

os.chdir(str(init_dir)+"/RESULTS/")

heptad_dir=(str(os.getcwd()))

for file in os.listdir(heptad_dir):
    #if file.endswith(".dat"):
    if str("THRESHOLD") not in file and str("TO_CLUSTER") not in file:
        #print("THIS IS FILE")
        #print(file)
        window_name=str(file).split("HEPTAD_ID_")
        #print("THIS IS FIRST NAME")
        first_part=str(window_name[1])
        #print(first_part)
        window_name2=str(first_part).split(".dat")
        #print("THIS IS SECOND NAME")
        second_part=str(window_name2[0])
       
        #window_name=str(file).split(".dat").str[0]
       
        #print("THIS IS COMPLETE NAME")
        current_df=pd.read_csv(file, sep="\s+")
        #print("THIS IS CURRENT DF")
        #print(current_df)
        column_names=current_df.columns
        #print("THESE ARE COLUMN NAMES")
        #print(column_names)
        current_df['b_score']=[-1.0 if v < 20 else 1.0 if v <=20 and v <30 else 1.5 if v <=30 and v<40 else 2.0 if v <=40 and v<50 else 3.0 if v <=50 and v<60 else 3.5 for v in current_df[column_names[2]]]
        current_df['c_score']=[-1.0 if v < 20 else 1.0 if v <=20 and v <30 else 1.5 if v <=30 and v<40 else 2.0 if v <=40 and v<50 else 3.0 if v <=50 and v<60 else 3.5 for v in current_df[column_names[3]]]
        current_df['e_score']=[-1.0 if v < 20 else 1.0 if v <=20 and v <30 else 1.5 if v <=30 and v<40 else 2.0 if v <=40 and v<50 else 3.0 if v <=50 and v<60 else 3.5 for v in current_df[column_names[5]]]
        current_df['f_score']=[-1.0 if v < 20 else 1.0 if v <=20 and v <30 else 1.5 if v <=30 and v<40 else 2.0 if v <=40 and v<50 else 3.0 if v <=50 and v<60 else 3.5 for v in current_df[column_names[6]]]
        current_df['g_score']=[-1.0 if v < 20 else 1.0 if v <=20 and v <30 else 1.5 if v <=30 and v<40 else 2.0 if v <=40 and v<50 else 3.0 if v <=50 and v<60 else 3.5 for v in current_df[column_names[7]]]
        current_df['STANDARD_HEPTAD_SCORE'] = current_df.apply(lambda row: row.b_score + row.c_score + row.e_score + row.f_score + row.g_score  , axis=1)
        current_df['ANTIGEN_NAME_ONLY']=current_df['emm_types'].str.split("regular").str[0]
        current_df.to_csv(str(cwd)+str("/HEPTAD_SITE_WITH_SCORE/")+str(second_part)+str("_all_heptad_scores.txt"), sep="\t")
        highest_heptad_score_for_each_antigen=current_df.sort_values('STANDARD_HEPTAD_SCORE').drop_duplicates('ANTIGEN_NAME_ONLY', keep='last')
        highest_heptad_score_for_each_antigen.to_csv(str(cwd)+str("/HEPTAD_SITE_WITH_SCORES_EACH ANTIGENS_BEST_HEPTAD_SCORING_WINDOW/")+str(second_part)+str("_best_windows_for_each_antigen.txt"), sep="\t") ### BEST HEPTAD SCORES FOR EACH ANTIGEN
        all_windows_with_score_at_least_10point5= current_df.loc[current_df['STANDARD_HEPTAD_SCORE'] >= 10.5]
        all_windows_with_score_at_least_10point5.to_csv(str(cwd)+str("/HEPTAD_SITE_WITH_SCORES_ALL_WINDOWS_GREATER_THAN_10.5/")+str(second_part)+str("_all_windows_scoring_at_least_10.5_.txt"), sep="\t")
    else:
        pass
os.chdir(str(init_dir)+"/PERCENTAGE_IDENTITY/")
percent_dir=(str(os.getcwd()))
for file in os.listdir(percent_dir):
    current_df=pd.read_csv(file, sep="\s+")
    percent_id_df=current_df.merge(current_df.Window.str.extractall(('(?<=IDENTICAL_)(.*?)(?=_regular)|(?<=>)(.*?)(?=_regular)')).reset_index(level=-1, drop=True), how='outer', left_index=True, right_index=True)       
    percent_id_df.columns= ['Window', 'Percent_Identity', 'first_match', 'second_match'] 
    percent_id_df['Antigen_Name_Only']= percent_id_df['first_match'].combine_first(percent_id_df['second_match'])
    del percent_id_df['first_match']
    del percent_id_df['second_match']
    percent_id_df=percent_id_df.sort_values('Percent_Identity').drop_duplicates('Antigen_Name_Only', keep='last')
    current_file_name=os.path.basename(file)
    current_df.to_csv(init_dir+str("/PERCENT_IDENTITY_PANDAS/")+str("ANTIGEN_NAME_ONLY_")+file, sep="\t")