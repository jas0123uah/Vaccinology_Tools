import os, glob, sys, os.path, copy, linecache, subprocess, re, collections
from statistics import mean 
import numpy as np


####################################  PLACE THINGS USER MUST EDIT HERE ####################################  

#fasta_file=open('/Users/jayspencer96/Downloads/test_2/new_fasta_file.fasta', 'r')
fasta_file=open('/Volumes/Mac_air/Mac_Air/M_prots_at_90/hole.txt', 'r')
#USE GETCWD
#register_file=open('/Users/jayspencer96/Downloads/test_2/new_register.txt', 'r')
register_file=open('/Volumes/Mac_air/Mac_Air/M_prots_at_90/reg.txt', 'r')

number_of_windows_to_get=96 ##total number of windows to use in your coiled coil proteins
number_of_proteins_to_build=4
number_of_windows_for_a_protein=24
PATH_TO_MARCOIL=""
MAIN_DIRECTORY=("/Volumes/Mac_air/Mac_Air/M_prots_at_90/")
PATH_TO_FINAL_OUTPUT=("/Users/jayspencer96/Desktop/90/")
window_delimiter=str(".") # THE DELIMITER WHICH CAN BE USED TO ENSURE EACH ANTIGEN IS ONLY COUNTED ONCE


####################################  PLACE THINGS USER MUST EDIT ABOVE ####################################
antigens_similar_to_vaccine=open(str(PATH_TO_FINAL_OUTPUT)+str("ANTIGEN_WINDOWS_WITH_HIGH_ID_TO_VACCINE_WINDOWS.txt"), "a+")
antigens_similar_to_vaccine.write(str("VACCINE WINDOW")+str("\t")+str("SIMILAR WINDOW") +str("\t")+str("Percent Identity")+str("\n"))

domain_length_pattern =re.compile(r'length = \d+\)')
d_number_pattern = re.compile(r'length = (.*?)\)')
percentage_pattern= re.compile(r"\d+\.\d")
header_pattern=re.compile(r'[^\s]+')
match_dict={}
list_of_windows_used=[]
drop_set=set()
count_for_suggested_proteins=0 #AN ITERATOR THAT WILL GO UP BY ONE WITH EACH SUGGESTED CC PROT APPENDED TO THE SUGGESTED PROTEINS FILE. 
full_path_best_windows=[]
basename_best_windows=[]

files = glob.glob(os.path.join(str(MAIN_DIRECTORY)+str("HIGH_PERCENTAGE_IDENTITY/"), "*.out"))
#files = glob.glob(os.path.join(str(MAIN_DIRECTORY)+str("HIGH_PERCENTAGE_IDENTITY/"), "*.out"))
###THIS WILL HAVE TO BE A SYS ARG /Users/jayspencer96/Downloads/test_2/ JUST USE GETCWD 

#USE GETCWD
heptad_list=["1","2","3","4", "5", "6", "7"]


#dict_with_full_file_path_and_actual_AA_seq_for_best_matches={}



list_of_headers=[]
for line in fasta_file:
    if ">" in line:
        line=line.strip()
        line=line[1:]
        line_to_add=str("high_id_w_")+str(line)+str(".out")
        list_of_headers.append(line_to_add)
        
    else:
        pass

print("THIS IS LIST OF HEADERS")
print(list_of_headers)
register_list=[]
sub_list=[]
index_iterator=0
for line in register_file:
    htad_letter=str(line.strip())
    path_to_aminos=list_of_headers[index_iterator]
    path_to_aminos=(str(path_to_aminos[10:-3] +str("faa")))
    path_to_aminos=str(MAIN_DIRECTORY)+ str(path_to_aminos)
    #path_to_aminos=str(MAIN_DIRECTORY)+str("/HEPTAD_SITE_FASTA/")+ str(path_to_aminos)
    current_aa_seq = linecache.getline(path_to_aminos, 2)
    current_aa_seq=current_aa_seq.strip()
    sub_list.append(current_aa_seq)
    sub_list.append(htad_letter)
    register_list.append(sub_list)
    sub_list=[]
    index_iterator+=1

print("THESE ARE THE HEPTAD REGISTERS")
print(register_list)

print("THIS IS THE WINDOW SIZE OF YOUR FRAGMENTS")

window_size=len(current_aa_seq)
print(window_size)
    #register_list.append(line.strip())
    

dict_of_windows_and_registers= dict(zip(list_of_headers, register_list))


#print("THIS IS DICT OF WDOWS")
#print(dict_of_windows_and_registers)











iterator_for_number_of_windows_to_get=0
while iterator_for_number_of_windows_to_get < number_of_windows_to_get: #THE NUMBER OF WINDOWS TO BUILD YOUR CC PROTEINS OUT OF
    match_dict={} ###CLEAR THIS OUT BECAUSE THE NUMBER OF MATCHES FOR EACH WINDOW SHOULD BE FEWER AND FEWER AS THE DROP SET EXPANDS
    for f in files: #ITERATE OVER EACH HIGH ID FILE IN THE DIRECTORY
        present_working_file=open(f,"r")
        list_of_matches_for_current_file=[]
        for line in present_working_file: # FOR EACH LINE IN THE HIGH ID FILE
            header_for_current_line=re.search(header_pattern, line) #GET THE NAME OF THE HIGH ID MATCH
            header_for_current_line=str(header_for_current_line.group(0))
            antigen_name_only=str((header_for_current_line.split(str(window_delimiter),1)[0])) #YOU DON'T WANT TO SAY A WINDOW IS COVERING MORE ANTIGENS THAN IT IS. FOR EXAMPLE, IF A PARTICULAR WINDOW HAS HIGH HEPTAD IDENTITY TO M1 WINDOW 1, IT LIKELKY HAS HIGH HEPTAD IDENTITY TO M1 WINDOW 2. M1 IS STILL ONLY ONE ANTIGEN SO YOU ONLY WANT TO COUNT IT ONCE. 
            if header_for_current_line != None and antigen_name_only not in drop_set and antigen_name_only not in list_of_matches_for_current_file: #MAKE SURE ITS NOT IN A DROP LIST FOR STUFF THATS ALREADY BEEN ACCOUNTED FOR BY ANOTHER WINDOW
                #print("THIS IS F")
                #print(f)
                #print("THIS IS HEADER FOR CURRENT LINE")
                #print(header_for_current_line)
                #print("THIS IS THE CURRENT FILE")
                #print(str(f))
                list_of_matches_for_current_file.append(antigen_name_only)
                #print("THESE ARE THE MATCHES FOR THE CURRENT FILE")
                #print(list_of_matches_for_current_file)

                #print("THIS IS ANTIGEN NAME ONLY")
                #print(antigen_name_only)
                #print("THIS THE DROP")
                #print(drop_set)
        #print("THIS IS F")
        #print(f)
        #print("THESE ARE F'S MATCHES")
        #print(list_of_matches_for_current_file)
        number_of_matches_after_drops=len(list_of_matches_for_current_file) #THE NUMBER OF ANTIGENS MATCHING HIGH ID FILE F.
        match_dict[str(f)] =number_of_matches_after_drops # A DICTIONARY WITH THE NAME OF THE FILE AND THE NUMBER OF ANTIGENS IT HAS HIGH HEPTAD ID WITH
        #print("THIS IS THE MATCH DICT")
        #print(match_dict)
    #test_job=open("/Users/jayspencer96/Desktop/match_dict.txt", "w")    
    #test_job.write(str(match_dict))
    #test_job.close()
    files_with_the_most_matches=[k for k, v in match_dict.items() if v == max(match_dict.values())] #https://stackoverflow.com/questions/47861737/how-to-get-multiple-max-key-values-in-a-dictionary
    #print("THESE ARE THE FILES WITH MOST MATCHES")
    #print(files_with_the_most_matches)
    #print("THIS IS THE DROP LIST")
    #print(drop_set)
    #input("PAUSE")
    #print("DUMBASS")
    #print(match_dict)
   # print("THE FILES THAT MATCH THE MOST ANTIGENS ARE:")
    #print(files_with_the_most_matches)
    
    #iterator_for_number_of_windows_to_get+=1
    #print("THE MATCH DICT IS")
   # print(match_dict)







    if len(files_with_the_most_matches)>1: #YOU MIGHT HAVE TWO OR MORE WINDOWS WITH THE MOST NUMBER OF MATCHES. GET THE BEST WINDOW BY LOOKING AT THE ACTUAL HEPTAD ID SCORES FOR THESE WINDOWS. HANDLING THE SAME ANTIGEN, DIFF WINDOWS GET HIGHEST SCORE?
        print("THE FOLLOWING FILES ALL HAVE THE SAME NUMBER OF MATCHES. LOOKING TO SEE WHICH OF THESE FILES HAS THE HIGHEST HEPTAD ID ")
        #print(files_with_the_most_matches)
        #input("WE NEED TO PAUSE")
        dictionary_for_files_with_most_matches={} #FILES WITH THE MOST MATCHES AND THEIR AVERAGE HEPTAD IDENTITY TO THE ANTIGENS THEY MATCH
        for file_with_lots_of_matches in files_with_the_most_matches:
            #list_of_percentages_for_current_file=[] SHOULD BE DICTIONARY OF PERCENTAGES BC YOU MIGHT HAVE A MATCH TO MORE THAN ONE WINDOW FROM THE SAME ANTIGEN. WANT THE HIGHEST HEPTAD ID FOR EACH ANTIGEN.
            dict_of_percentages_for_current_file={} #GET THE HIGHEST HEPTAD IDENTITY FOR EACH MATCHED ANTIGEN IN A FILE. THIS DICTIONARY WILL BE USED TO CALCULATE THE AVERAGE HEPTAD IDENTITY FOR A PARTICULAR FILE/WINDOW.
            current_file=open(file_with_lots_of_matches, "r")
            for line in current_file: #GET PERCENT HEPTAD ID FOR EACH MATCH
                percentage_for_current_line=re.search(percentage_pattern, line)
                percentage_for_current_line=float(percentage_for_current_line.group(0))
                header_for_current_line=re.search(header_pattern, line) #GET THE NAME OF THE HIGH ID MATCH
                header_for_current_line=str(header_for_current_line.group(0))
                antigen_name_only=str((header_for_current_line.split(str(window_delimiter),1)[0]))
                if antigen_name_only not in dict_of_percentages_for_current_file:
                    dict_of_percentages_for_current_file[str(antigen_name_only)]=float(percentage_for_current_line)
                elif antigen_name_only in dict_of_percentages_for_current_file:
                    current_highest_percentage=dict_of_percentages_for_current_file.get(antigen_name_only)
                    if float(current_highest_percentage) >= float(percentage_for_current_line):
                        pass
                    elif float(current_highest_percentage) < float(percentage_for_current_line):
                        dict_of_percentages_for_current_file[antigen_name_only]= float(percentage_for_current_line)

            
            average_for_current_file=np.array(list(dict_of_percentages_for_current_file.values())).mean()

            #average_for_current_file=mean(dict_of_percentages_for_current_file) #FIND THE AVERAGE HEPTAD ID OF THE MATCHES
            dictionary_for_files_with_most_matches[str(file_with_lots_of_matches)]=average_for_current_file # A DICTIONARY WITH THE FILE NAMES AND AVERAGE HEPTAD ID


        most_matches_and_highest_avg=[k for k, v in dictionary_for_files_with_most_matches.items() if v == max(dictionary_for_files_with_most_matches.values())]
        print("THIS IS MOST MATCHES AND HIGHEST AVG")
        print(most_matches_and_highest_avg)
        #print("THE FILE WITH THE HIGHEST HEPTAD ID AVG FROM THE DICT WITH THE MOST MATCHES IS: ")
        #print(most_matches_and_highest_avg[0]) # AT THIS POINT IF THE LENGTH IS GREATER THAN ONE WE'VE LOOKED FAR ENOUGH, JUST GET THE FIRST WINDOW IF THERE IS STILL MORE THAN ONE
        full_path_best_windows.append(most_matches_and_highest_avg[0])
        basename=os.path.basename(most_matches_and_highest_avg[0])
        basename_best_windows.append(basename)
        current_best_window= str(most_matches_and_highest_avg[0])
    
    elif len(files_with_the_most_matches)==1:
        print("I AM HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
        full_path_best_windows.append(files_with_the_most_matches[0])
        basename=os.path.basename(files_with_the_most_matches[0])
        basename_best_windows.append(basename)
        current_best_window= str(files_with_the_most_matches[0])
        




    most_matched_file=open(current_best_window, "r")
    print("THIS IS MOST MATCHED FILE")
    print(str(current_best_window))
    


    for line in most_matched_file: #WE NEED TO KNOW WHAT WE HAVE COVERED IN OUR FINAL SELECTION. ADD THESE WINDOWS TO A DROP LIST
        #print("THIS IS LINE")
        #print(line)
        header=re.search(header_pattern, line)
        header=str(header.group(0))
        #print("THIS IS HEADER")
        #print(header)
        if header != None:
            antigen_name_only=str((header.split(str(window_delimiter),1)[0])) #YOU DON'T WANT TO SAY A WINDOW IS COVERING MORE ANTIGENS THAN IT IS. FOR EXAMPLE, IF A PARTICULAR WINDOW HAS HIGH HEPTAD IDENTITY TO M1 WINDOW 1, IT LIKELKY HAS HIGH HEPTAD IDENTITY TO M1 WINDOW 2. M1 IS STILL ONLY ONE ANTIGEN SO YOU ONLY WANT TO COUNT IT ONCE. 
            #print("THIS IS ANTIGEN NAME ONLY")
            #print(antigen_name_only)
            drop_set.add(antigen_name_only)
           
            everything_before_HIGH_PERCENTAGE_IDENTITY=current_best_window.split('HIGH_PERCENTAGE_IDENTITY')[0]
           
            file_path_to_actual_AA_seq=str(everything_before_HIGH_PERCENTAGE_IDENTITY)+str("HEPTAD_SITE_FASTA/") +str(basename[10:-3] +str("faa"))
            #print("THIS IS FILE PATH TO AA SEQ")
            #print(file_path_to_actual_AA_seq)
            
            basename_for_file_containing_AA_seq=str(basename[10:-3] +str("faa"))
            
                
            #basename_for_file_containing_AA_seq=str(basename[10:-3] +str("faa"))
            
            
        else:
            pass
    #print("THIS IS THE CURRENT DROP SET")
    #print(drop_set)
    iterator_for_number_of_windows_to_get +=1  #### AT THE END OF THIS WHILE LOOP YOU SHOULD HAVE YOUR MOST CONSERVED HEPTAD WINDOWS

#print("THIS IS DROP SET")
#print(drop_set)
#print("THIS IS FULL PATH BEST WINDOWS")
#print(full_path_best_windows)
dict_of_best_windows_and_their_registers = {key: dict_of_windows_and_registers[key] for key in dict_of_windows_and_registers.keys() & basename_best_windows} #CREATE A SMALL DICTIONARY WITH YOUR BEST WINDOWS AND THEIR HEPTAD LETTER. RIGHT NOW THIS IS JUST THE BEST WINDOW FASTA HEADERS AND THEIR LETTER AS VALUE. NEED TO GET THE SEQUENCE ASSOCIATED WITH THE HEADER

#print("THIS IS BASENAME BEST WINDOWS")
#print(basename_best_windows)

#print("THIS IS DICT OF WINDOWS AND REGISTERS")
#print(dict_of_windows_and_registers)
print("THIS IS DICT OF BEST WINDOWS AND REGISTERS")
print(dict_of_best_windows_and_their_registers)








#for k, v in dict_of_best_windows_and_their_registers.items():
 #   print(k)
  #  print(v[1])

a_windows=[]
b_windows=[]
c_windows=[]
d_windows=[]
e_windows=[]
f_windows=[]
g_windows=[]

#print("THIS IS B DICT")
for item in dict_of_best_windows_and_their_registers.items():
    #print("THIS IS ITEM 1.1")
    #print(item[1][1])
    if item[1][1] == str("1"):
        a_windows.append(item)
    elif item[1][1] ==str("2"):
        b_windows.append(item)
    elif item[1][1] ==str("3"):
        c_windows.append(item)
    elif item[1][1] ==str("4"):
        d_windows.append(item)
    elif item[1][1] ==str("5"):
        e_windows.append(item)
    elif item[1][1] ==str("6"):
        f_windows.append(item)
    elif item[1][1] ==str("7"):
        g_windows.append(item)


stored_a_windows=copy.deepcopy(a_windows)
stored_b_windows=copy.deepcopy(b_windows)
stored_c_windows=copy.deepcopy(c_windows)
stored_d_windows=copy.deepcopy(d_windows)
stored_e_windows=copy.deepcopy(e_windows)
stored_f_windows=copy.deepcopy(f_windows)
stored_g_windows=copy.deepcopy(g_windows)    



print("THIS IS A WINDOWS")
print(a_windows)
print("THIS IS B WINDOWS")
print(b_windows)
print("THIS IS C WINDOWS")
print(c_windows)
print("THIS IS D WINDOWS")
print(d_windows)
print("THIS IS E WINDOWS")
print(e_windows)
print("THIS IS F WINDOWS")
print(f_windows)
print("THIS IS G WINDOWS")
print(g_windows)


################INSERT CODE HERE TO WRITE THE NAMES OF THE A,B,C,D,E,F, AND G WINDOWS TO A SELECTED A WINDOW, B WINDOW FILE





remainder=(window_size%7) #CALCULATE HOW MANY COMPLETE A-G CYCLES EACH WINDOW MAKES AND WHICH LETTERS ARE AFTER. MINUS 1 IS USED TO GET THE NUMBER OF HEPTAD LETTERS TO ADVANCE TO GET THE FINAL HEPTAD LETTER WHICH FRAGS GET PUT TOGETHER DEPENDS ON THE REMAINDER SIZE!!!!

print("THIS IS REMAINDER")
print(remainder)
ideal_CC_prot_version_A=[]
windows_used_for_ideal_CC_prot_version_A=[]



























initial_letter=1
if remainder ==1:
    
    last_window_used=initial_letter
    number_of_residues_trimmed=0
    iterator=0
    dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values= collections.OrderedDict()
#######ADD WHILE LOOP HERE TO BUILD N NUMBER OF PROTS. DROP FROM STORED_WINDOWS.
    number_of_proteins_made=0
    while number_of_proteins_made < number_of_proteins_to_build:
        
        
    
        for key in list_of_windows_used: #REMOVES THE WINDOWS WHICH WERE USED TO MAKE THE PREVIOUS PROTEIN
           
            for second_key in stored_a_windows:
                if key==second_key[0]:
                    stored_a_windows.remove(second_key)
            for second_key in stored_b_windows:
                if second_key[0]== key:
                    stored_b_windows.remove(second_key)
            for second_key in stored_c_windows:
                if second_key[0]== key:
                    stored_c_windows.remove(second_key)
            for second_key in stored_d_windows:
                if second_key[0]== key:
                    stored_d_windows.remove(second_key)
            for second_key in stored_e_windows:
                if second_key[0]== key:
                    stored_e_windows.remove(second_key)
            for second_key in stored_f_windows:
                if second_key[0]== key:
                    stored_f_windows.remove(second_key)
            for second_key in stored_g_windows:
                if second_key[0]== key:
                    stored_g_windows.remove(second_key)
               
        print("THESE ARE THE WINDOWS REMAINING AFTER THE PREVIOUS ITERATION")
        print("REMAINING A WINDOWS")
        print(stored_a_windows)
        print("REMAINING B WINDOWS")
        print(stored_b_windows)
        print("REMAINING C WINDOWS")
        print(stored_c_windows)
        print("REMAINING D WINDOWS")
        print(stored_d_windows)
        print("REMAINING E WINDOWS")
        print(stored_e_windows)
        print("REMAINING F WINDOWS")
        print(stored_f_windows)
        print("REMAINING G WINDOWS")
        print(stored_g_windows)
        
                
          
            
        while initial_letter <= 7: #USE THIS WHILE LOOP TO BUILD ALL 7 VARIANTS OF YOUR NTH IDEAL CC PROT 
            last_window_used = initial_letter
            a_windows=copy.deepcopy(stored_a_windows) #WINDOWS GET POPPED WITH EACH ITERATION TO BUILD THE POSSIBLE IDEAL CC. YOU NEED TO BE ABLE TO RESTORE THE LIST BY CALLING ON THE STORED LIST
            b_windows=copy.deepcopy(stored_b_windows)
            c_windows=copy.deepcopy(stored_c_windows)
            d_windows=copy.deepcopy(stored_d_windows)
            e_windows=copy.deepcopy(stored_e_windows)
            f_windows=copy.deepcopy(stored_f_windows)
            g_windows=copy.deepcopy(stored_g_windows)
            print("YEE HAW")
            print(iterator)

        

            while iterator < number_of_windows_for_a_protein: ####MAKE SYS ARG
                if last_window_used== 1 or last_window_used== str(1) : #THE VALUE FOR A START WITH A AND END WITH A 
                    try: ### See if B list windows are empty WE ARE TRYING TO BUILD A PROTEIN WITH B AS THE FIRST REGISTER. TRY B LIST WINDOWS YOU TRIM NO RESIDUES
                        ideal_CC_prot_version_A.append(b_windows[0][1][0])
                        last_window_used=heptad_list[0]
                        fasta_header=str(a_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        b_windows.pop(0)
                        iterator+=1
                    except IndexError: ###  YOU NEED TO START WITH AN "B" LETTER.  "B" list windows are empty, try A list windows. You only have to trim 1 residue 
                        try:
                            current_frag=a_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[6]
                            fasta_header=str(a_windows[0][0])
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            a_windows.pop(0)
                            iterator+=1  
                        except IndexError: ### YOU NEED TO START WITH A "B" LETTER. B and A list windows are empty, try G list windows. You trim 2 residues
                            try:
                                current_frag=g_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[5]
                                fasta_header=str(g_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                g_windows.pop(0)
                                iterator+=1
                            except IndexError: ### YOU NEED TO START WITH A "B" LETTER. B, A, AND G list windows are empty, try F list windows. You trim 3 residues
                                try:
                                    current_frag=f_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[4]
                                    fasta_header=str(f_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    f_windows.pop(0)
                                    iterator+=1
                                except IndexError: ###### B,A,G, AND F  list windows are empty, try E list windows. You trim 4 residues
                                    try:
                                        current_frag=e_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[3]
                                        fasta_header=str(e_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        e_windows.pop(0)
                                        iterator+=1
                                    except IndexError:
                                        try: ###### B,A,G,F,AND E list windows are empty, try D list windows You trim 5 residues
                                            current_frag=d_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag)
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[2]
                                            fasta_header=str(d_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            d_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ###### B,A,G,F,E, AND D list windows are empty, try C list windows You trim 6 residues
                                                current_frag=b_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[1]
                                                fasta_header=str(b_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                b_windows.pop(0)
                                                iterator+=1
                if last_window_used ==2 or last_window_used== str(2): ### WE JUST USED A "B" GET A "C" 
                    try:
                        ideal_CC_prot_version_A.append(c_windows[0][1][0])
                        last_window_used=heptad_list[1]
                        fasta_header=str(c_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        c_windows.pop(0)
                        iterator+=1
                    except IndexError: ### C list windows are empty, try B list windows You trim 1 residue
                        try:
                            current_frag=b_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[0]
                            fasta_header=str(b_windows[0][0])
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            b_windows.pop(0)
                            iterator+=1
                        except IndexError: ### C AND B list windows are empty, try A list windows You trim 2 residues
                            try:
                                current_frag=a_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[6]
                                fasta_header=str(a_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                a_windows.pop(0)
                                iterator+=1
                            except IndexError: ###### C,B,AND A list windows are empty, try G list windows You trim 3 residues
                                try:
                                    current_frag=g_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[5]
                                    fasta_header=str(g_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    g_windows.pop(0)
                                    iterator+=1
                                except IndexError:
                                    try: ###### C,B,A,AND G list windows are empty, try F list windows. You trim 4 residues
                                        current_frag=f_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[4]
                                        fasta_header=str(f_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        f_windows.pop(0)
                                        iterator+=1
                                    except IndexError:
                                        try: ###### C,B,A,G AND F list windows are empty, try E list windows. You trim 5 residues
                                            current_frag=e_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag)
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[3]
                                            fasta_header=str(e_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            e_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ### C,B,A,G,F, AND E lists are empty, try D list You trim 6 residues
                                                current_frag=d_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[2]
                                                fasta_header=str(d_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                d_windows.pop(0)
                                                iterator+=1
                if last_window_used ==3 or last_window_used== str(3):
                    try:
                        ideal_CC_prot_version_A.append(d_windows[0][1][0])
                        last_window_used=heptad_list[2]
                        fasta_header=str(d_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        d_windows.pop(0)
                        iterator+=1
                    except IndexError: ### D list windows are empty, try C list windows You trim 1 residue
                        try:
                            current_frag=c_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[1]
                            fasta_header=str(c_windows[0][0])
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            c_windows.pop(0)
                            iterator+=1
                        except IndexError: ### D AND C list windows are empty, try B list windows You trim 2 residues
                            try:
                                current_frag=b_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[0]
                                fasta_header=str(b_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                b_windows.pop(0)
                                iterator+=1
                            except IndexError: ###### D,C,AND B list windows are empty, try A list windows You trim 3 residues
                                try:
                                    current_frag=a_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[6]
                                    fasta_header=str(a_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    a_windows.pop(0)
                                    iterator+=1
                                except IndexError:
                                    try: ###### D,C,B,AND A list windows are empty, try G list windows. You trim 4 residues
                                        current_frag=g_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[5]
                                        fasta_header=str(g_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        g_windows.pop(0)
                                        iterator+=1 
                                    except IndexError:
                                        try: ###### D,C,B,A AND G list windows are empty, try F list windows. You trim 5 residues
                                            current_frag=f_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag)
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[4]
                                            fasta_header=str(f_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            f_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ### D,C,B,A,G, AND F lists are empty, try E list You trim 6 residues
                                                current_frag=e_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[3]
                                                fasta_header=str(e_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                e_windows.pop(0)
                                                iterator+=1


                if last_window_used ==4 or last_window_used== str(4): ### WE JUST USED A "D" GET AN "E"
                    try:
                        ideal_CC_prot_version_A.append(e_windows[0][1][0])
                        last_window_used=heptad_list[3]
                        fasta_header=str(e_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        e_windows.pop(0)
                        iterator+=1
                    except IndexError: ### E list windows are empty, try D list windows You trim 1 residue
                        try:
                            current_frag=d_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[2]
                            fasta_header=str(d_windows[0][0])
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            d_windows.pop(0)
                            iterator+=1
                        except IndexError: ### E AND D list windows are empty, try C list windows You trim 2 residues
                            try:
                                current_frag=c_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[1]
                                fasta_header=str(c_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                c_windows.pop(0)
                                iterator+=1
                            except IndexError: ###### E,D, AND C list windows are empty, try B list windows You trim 3 residues
                                try:
                                    current_frag=b_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[0]
                                    fasta_header=str(b_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    b_windows.pop(0)
                                    iterator+=1
                                except IndexError:
                                    try: ###### E,D,C,AND B list windows are empty, try A list windows. You trim 4 residues
                                        current_frag=a_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[6] 
                                        fasta_header=str(a_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        a_windows.pop(0)
                                        iterator+=1
                                    except IndexError:
                                        try: ###### E,D,C,B, AND A list windows are empty, try G list windows. You trim 5 residues
                                            current_frag=g_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag) 
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[5]
                                            fasta_header=str(g_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            g_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ### E,D,C,B,A AND G lists are empty, try F list You trim 6 residues
                                                current_frag=f_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[4]
                                                fasta_header=str(f_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                f_windows.pop(0)
                                                iterator+=1
                if last_window_used ==5 or last_window_used== str(5): ### WE JUST USED A "E" GET AN "F"
                    try:
                        ideal_CC_prot_version_A.append(f_windows[0][1][0])
                        last_window_used=heptad_list[4]
                        fasta_header=str(f_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        f_windows.pop(0)
                        iterator+=1
                    except IndexError: ### F list windows are empty, try E list windows You trim 1 residue
                        try:
                            current_frag=e_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[3]
                            fasta_header=str(e_windows[0][0]) 
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            e_windows.pop(0)
                            iterator+=1
                        except IndexError: ### F AND E list windows are empty, try D list windows You trim 2 residues
                            try:
                                current_frag=d_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[2]
                                fasta_header=str(d_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                d_windows.pop(0)
                                iterator+=1
                            except IndexError: ###### F,E,AND D list windows are empty, try C list windows You trim 3 residues
                                try:
                                    current_frag=c_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[1]
                                    fasta_header=str(c_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    c_windows.pop(0)
                                    iterator+=1
                                except IndexError:
                                    try: ###### F,E,D,AND C list windows are empty, try B list windows. You trim 4 residues
                                        current_frag=b_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[0] 
                                        fasta_header=str(b_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        b_windows.pop(0)
                                        iterator+=1
                                    except IndexError:
                                        try: ###### F,E,D,C, AND B list windows are empty, try A list windows. You trim 5 residues
                                            current_frag=a_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag)
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[6]
                                            fasta_header=str(a_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            a_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ### F,E,D,C,B, AND A lists are empty, try G list You trim 6 residues
                                                current_frag=g_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[5]
                                                fasta_header=str(g_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                g_windows.pop(0)
                                                iterator+=1
                if last_window_used ==6 or last_window_used== str(6): ### WE JUST USED A "F" GET A "G"
                    try:
                        ideal_CC_prot_version_A.append(g_windows[0][1][0])
                        last_window_used=heptad_list[5]
                        fasta_header=str(g_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        g_windows.pop(0)
                        iterator+=1
                    except IndexError: ### G list windows are empty, try F list windows You trim 1 residue
                        try:
                            current_frag=f_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[4]
                            fasta_header=str(f_windows[0][0])
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            f_windows.pop(0)
                            iterator+=1
                        except IndexError: ### G, AND F list windows are empty, try E list windows You trim 2 residues
                            try:
                                current_frag=d_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[3]
                                fasta_header=str(d_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                d_windows.pop(0)
                                iterator+=1
                            except IndexError: ###### F,E,AND D list windows are empty, try D list windows You trim 3 residues
                                try:
                                    current_frag=c_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[2]
                                    fasta_header=str(c_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    c_windows.pop(0)
                                    iterator+=1
                                except IndexError:
                                    try: ###### F,E,D,AND C list windows are empty, try C list windows. You trim 4 residues
                                        current_frag=b_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[1]
                                        fasta_header=str(b_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        b_windows.pop(0)
                                        iterator+=1
                                    except IndexError:
                                        try: ###### F,E,D,C, AND B list windows are empty, try B list windows. You trim 5 residues
                                            current_frag=a_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag)
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[0]
                                            fasta_header=str(a_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            a_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ### F,E,D,C,B, AND A lists are empty, try A list You trim 6 residues
                                                current_frag=g_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[6]
                                                fasta_header=str(g_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                g_windows.pop(0)
                                                iterator+=1


                if last_window_used == 7 or last_window_used== str(7): ### WE JUST USED A "G" GET ANOTHER "G"
                    try:
                        ideal_CC_prot_version_A.append(g_windows[0][1][0])
                        last_window_used=heptad_list[6]
                        fasta_header=str(g_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        g_windows.pop(0)
                        iterator+=1
                    except IndexError: ### G list windows are empty, try F list windows You trim 1 residue
                        try:
                            current_frag=f_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[5]
                            fasta_header=str(f_windows[0][0])
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            f_windows.pop(0)
                            iterator+=1
                        except IndexError: ### G AND F list windows are empty, try E list windows You trim 2 residues
                            try:
                                current_frag=e_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[4]
                                fasta_header=str(e_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                e_windows.pop(0)
                                iterator+=1
                            except IndexError: ###### G,F,AND E list windows are empty, try D list windows You trim 3 residues
                                try:
                                    current_frag=d_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[3]
                                    fasta_header=str(d_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    d_windows.pop(0)
                                    iterator+=1
                                except IndexError:
                                    try: ###### G,F,E,AND D list windows are empty, try C list windows. You trim 4 residues
                                        current_frag=c_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[2]
                                        fasta_header=str(c_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        c_windows.pop(0)
                                        iterator+=1
                                    except IndexError:
                                        try: ###### G,F,E,D AND C list windows are empty, try B list windows. You trim 5 residues
                                            current_frag=b_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag)
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[1]
                                            fasta_header=str(b_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            b_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ### G,F,E,D,C, AND B lists are empty, try A list You trim 6 residues
                                                current_frag=c_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[0]
                                                fasta_header=str(a_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                a_windows.pop(0)
                                                iterator+=1















            #nYAYYY
            final_ideal_CC_prot="".join(ideal_CC_prot_version_A)
            windows_used_for_ideal_CC_prot_version_A.append(float(number_of_residues_trimmed))
            print("THIS THE FINAL PROT")
            print(final_ideal_CC_prot)
            #input("PAUSE")
            dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values[str(final_ideal_CC_prot)]=windows_used_for_ideal_CC_prot_version_A
            ideal_CC_prot_version_A=[] #RESET TO BUILD CC PROT STARTING WITH WINDOW B, THEN WINDOW C, THEN WINDOW D...
            number_of_residues_trimmed=0 #RESET THE TRIM RESIDUE COUNT
            iterator=0 #RESET TO 0 SO YOU CAN GET THE K HEPTAD WINDOWS, FOR YOUR NTH CC PROTEIN 
            windows_used_for_ideal_CC_prot_version_A=[] #RESET TO GET WINDOWS USED TO BUILD CC PROT STARTING WITH WINDOW B, THEN WINDOW C, THEN WINDOW D...

            initial_letter+=1
        print("THIS IS THE DICT OF CC PROTS CREATED")
        print(dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values)
        #CC_prot_with_fewest_residues_trimmed = min(dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values.values()) 

        min_value = min(dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values.values(), key=lambda x: x[-1])[-1]
        dict_of_CC_proteins_w_fewest_residues_trimmed = {k: v for k, v in dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values.items() if v[-1] == min_value}
        print("THESE ARE THE PROTEINS WITH FEWEST RESIDUES TRIMMED")
        print(dict_of_CC_proteins_w_fewest_residues_trimmed)
        
        if len(dict_of_CC_proteins_w_fewest_residues_trimmed.keys()) >1:
            
            print("THERE IS MORE THAN ONE CC PROTEIN WHICH HAS THE FEWEST RESIDUES TRIMMED. WRITING THESE PROTEINS TO A TEMPORARY TEXT FILE AND RUNNING MARCOIL ANALYSIS TO SELECT THE PROTEIN WITH THE HIGHEST COILED COIL PROBABILITY.")
            temp_cc_prot_iterator=0 #SET TO 0 SO YOU CAN USE THE NUMBER IN THE FASTA HEADER TO EXTRACT BY KEY INDEX
            print(os.getcwd())
            os.chdir("/Users/jayspencer96/MARCOIL") #####MAKE SYS ARG
            print(os.getcwd())
            temporary_text_file_for_MARCOIL_analysis=open('temp.txt', 'w')
            for key in dict_of_CC_proteins_w_fewest_residues_trimmed:
                print("THIS IS KEY")
                print(key)
                temporary_text_file_for_MARCOIL_analysis.write(">CC_protein_" +str(temp_cc_prot_iterator)+str("\n"))
                temporary_text_file_for_MARCOIL_analysis.write(str(key) +str("\n"))
                temp_cc_prot_iterator+=1
            temporary_text_file_for_MARCOIL_analysis.close()

            os.system("./marcoil +dlsS temp.txt")
            os.chdir('/Users/jayspencer96/MARCOIL/Outputs/')###NOTE output_for_MARCOIL must be an alias in your bash_profile that navigates to MARCOIL output
            print(os.getcwd())
            domains_file=open('Domains', 'r')




            dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values={}
            while True:

                        for line in domains_file:
                            if ">" in line:
                                current_CC_protein=str(line.strip())
                                list_of_domain_lengths=[]

                            elif "NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 2.0 :" in line:
                                number_of_2_percent_domains= int(line[-2]) ###Use this to determine how many lines you have to iterate through
                                if number_of_2_percent_domains == 0: ### No 2 percent domains. 

                                    list_of_domain_lengths.append(0)

                                elif number_of_2_percent_domains >=1:

                                    lines_to_iterate_through=1
                                    total_2_percent_domain_length=0
                                    while lines_to_iterate_through <= number_of_2_percent_domains:
                                        #print("THIS is LINE")
                                        #print(line)
                                        line=next(domains_file)
                                        #print("THIS IS NEXT LINE")
                                        #print(line)
                                        current_2_percent_domain_length=str(re.search(domain_length_pattern, line).group(0)) #GET THE NAME OF THE HIGH ID MATCH
                                        current_2_percent_domain_length = int(re.search(d_number_pattern, current_2_percent_domain_length).group(1))
                                        total_2_percent_domain_length=total_2_percent_domain_length+current_2_percent_domain_length
                                        lines_to_iterate_through+=1
                                    list_of_domain_lengths.append(total_2_percent_domain_length)



                            elif "NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 10.0 :" in line:
                                number_of_10_percent_domains= int(line[-2]) ###Use this to determine how many lines you have to iterate through

                                if number_of_10_percent_domains == 0: ### No 10 percent domains. 
                                    list_of_domain_lengths.append(0)


                                elif number_of_10_percent_domains >=1:
                                    lines_to_iterate_through=1
                                    total_10_percent_domain_length=0
                                    while lines_to_iterate_through <= number_of_10_percent_domains:
                                        line=next(domains_file)
                                        current_10_percent_domain_length=str(re.search(domain_length_pattern, line).group(0)) #GET THE NAME OF THE HIGH ID MATCH
                                        current_10_percent_domain_length = int(re.search(d_number_pattern, current_10_percent_domain_length).group(1))
                                        #current_10_percent_domain_length=int(current_10_percent_domain_length.group(0))
                                        total_10_percent_domain_length=total_10_percent_domain_length+current_10_percent_domain_length
                                        lines_to_iterate_through+=1
                                    list_of_domain_lengths.append(total_10_percent_domain_length)
                            elif "NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 50.0 :" in line:
                                number_of_50_percent_domains= int(line[-2]) ###Use this to determine how many lines you have to iterate through

                                if number_of_50_percent_domains == 0: ### No 50 percent domains. Go to the next line

                                    list_of_domain_lengths.append(0)


                                elif number_of_50_percent_domains >=1:
                                    lines_to_iterate_through=1
                                    total_50_percent_domain_length=0
                                    while lines_to_iterate_through <= number_of_50_percent_domains:
                                        line=next(domains_file)
                                        current_50_percent_domain_length=str(re.search(domain_length_pattern, line).group(0))  #GET THE NAME OF THE HIGH ID MATCH
                                        current_50_percent_domain_length = int(re.search(d_number_pattern, current_50_percent_domain_length).group(1))
                                        total_50_percent_domain_length=total_50_percent_domain_length+current_50_percent_domain_length
                                        lines_to_iterate_through+=1
                                    list_of_domain_lengths.append(total_50_percent_domain_length)
                            elif "NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 90.0 :" in line:
                                number_of_90_percent_domains= int(line[-2]) ###Use this to determine how many lines you have to iterate through

                                if number_of_90_percent_domains == 0: ### No 90 percent domains. Go to the next line
                                    list_of_domain_lengths.append(0)


                                elif number_of_90_percent_domains >=1:
                                    lines_to_iterate_through=1
                                    total_90_percent_domain_length=0
                                    while lines_to_iterate_through <= number_of_90_percent_domains:
                                        line=next(domains_file)
                                        current_90_percent_domain_length=str(re.search(domain_length_pattern, line).group(0))  #GET THE NAME OF THE HIGH ID MATCH
                                        current_90_percent_domain_length = int(re.search(d_number_pattern, current_90_percent_domain_length).group(1))
                                        total_90_percent_domain_length=total_90_percent_domain_length+current_90_percent_domain_length
                                        lines_to_iterate_through+=1
                                    list_of_domain_lengths.append(total_90_percent_domain_length)

                            elif "NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 99.0 :" in line:
                                number_of_99_percent_domains= int(line[-2]) ###Use this to determine how many lines you have to iterate through

                                if number_of_99_percent_domains == 0: ### No 10 percent domains. Go to the next line
                                    list_of_domain_lengths.append(0)
                                    dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values[str(current_CC_protein)] =list_of_domain_lengths

                                elif number_of_99_percent_domains >=1:
                                    lines_to_iterate_through=1
                                    total_99_percent_domain_length=0
                                    while lines_to_iterate_through <= number_of_99_percent_domains:
                                        line=next(domains_file)
                                        current_99_percent_domain_length=str(re.search(domain_length_pattern, line).group(0))  #GET THE NAME OF THE HIGH ID MATCH
                                        current_99_percent_domain_length = int(re.search(d_number_pattern, current_99_percent_domain_length).group(1))
                                        total_99_percent_domain_length=total_99_percent_domain_length+current_99_percent_domain_length
                                        lines_to_iterate_through+=1
                                    list_of_domain_lengths.append(total_99_percent_domain_length)
                                    dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values[str(current_CC_protein)] =list_of_domain_lengths

                            else:
                                pass
                        print("THIS IS VACCINE PROTS AND DOMAINS")
                        print(dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values)
                        break


            avg_dict = {k: (sum(v) / len(v)) for k, v in dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values.items()}

            print(avg_dict)
            best_protein=max(avg_dict, key=avg_dict.get)
            key_for_best_protein=int(int(best_protein[-1])-1)

            print("THIS IS AVG DICT")
            print(avg_dict)

            print("THIS PROTEIN HAS THE HIGHEST CC AVG SO IT WILL BE USED") #THERE MIGHT BE MORE THAN ONE, BUT AT THIS POINT YOU'VE COMPARED ENOUGH
            print(max(avg_dict, key=avg_dict.get))
            print("THIS IS THE INDEX NUMBER FOR THE BEST PROTEIN")
            print(key_for_best_protein)
            #items = list(ordered.items())
            #items = list(dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values.items())
            items = list(dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values.items())
            #FINAL_BEST_DICT=(items[key_for_best_protein])
            #print(items[dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values[key_for_best_protein]])
            final_dict={} ### START FRESH SO THERE IS ONLY EVER 1 KEY IN THIS DICT, AND THAT'S THE BEST PROTEIN
            final_dict=items[key_for_best_protein]  #THIS WORKSSSSSS

            print("THIS IS THE BEST IDEALIZED COILED COIL PROTEIN")
            print(final_dict[0])
            suggested_proteins_file=open(PATH_TO_FINAL_OUTPUT+str("Suggested_proteins.txt"), "a+")
            print(len(final_dict[1]))
            i=0

            suggested_proteins_file.write(">Suggested_Coiled_Coil_Protein_"+str(count_for_suggested_proteins)+str("\n"))
            suggested_proteins_file.write(str(final_dict[0])+str("\n"))
            count_for_suggested_proteins+=1
            suggested_proteins_file.close()
            windows_for_suggested_protein=open(PATH_TO_FINAL_OUTPUT+str("Windows_used_to_make_suggested_proteins.txt"),"a+")
            windows_for_suggested_protein.write("********** THESE ARE THE WINDOWS USED TO MAKE SUGGESTED PROTEIN "+str(count_for_suggested_proteins) + str(". ")+str(int((final_dict[1][-1])))+ str(" RESIDUES HAD TO BE TRIMMED IN TOTAL TO CREATE THE FINAL PROTEIN. ********** \n" ))
            print("THESE ARE THE WINDOWS USED TO BUILD YOUR IDEALIZED COILED COIL PROTEIN:")
            while i< len(final_dict[1]) - 1:
                window=str(final_dict[1][i])
                print(window[10:-4])
                windows_for_suggested_protein.write(">"+str(window[10:-4])+str("\n"))
                i+=1
            windows_for_suggested_protein.close()
            list_of_windows_used=final_dict[1][:-1]
            #print("THIS IS LIST OF WINDOWS USED")
            print(list_of_windows_used)
        
            number_of_proteins_made+=1
        
        
        
        
        
        
        elif len(dict_of_CC_proteins_w_fewest_residues_trimmed.keys()) == 1:
            print("THERE IS ONLY ONE CC PROTEIN WHICH HAS THE FEWEST RESIDUES TRIMMED.")
            final_dict={} ### START FRESH SO THERE IS ONLY EVER 1 KEY IN THIS DICT, AND THAT'S THE BEST PROTEIN
            #final_dict=items[key_for_best_protein]  #THIS WORKSSSSSS

            print("THIS IS THE BEST IDEALIZED COILED COIL PROTEIN")
            my=list(dict_of_CC_proteins_w_fewest_residues_trimmed.items())
            print(my[0][0]) #THE LIST OF WINDOWS USED AND THE NUM RES TRIMMED AT END
            #print(my[0][1][i])
            
            suggested_proteins_file=open(PATH_TO_FINAL_OUTPUT+str("Suggested_proteins.txt"), "a+")
            
            i=0

            suggested_proteins_file.write(">Suggested_Coiled_Coil_Protein_"+str(count_for_suggested_proteins)+str("\n"))
            suggested_proteins_file.write(str(my[0][0])+str("\n"))
            count_for_suggested_proteins+=1
            suggested_proteins_file.close()
            windows_for_suggested_protein=open(PATH_TO_FINAL_OUTPUT+str("Windows_used_to_make_suggested_proteins.txt"),"a+")
            windows_for_suggested_protein.write("********** THESE ARE THE WINDOWS USED TO MAKE SUGGESTED PROTEIN "+str(count_for_suggested_proteins) + str(". ")+str(int((my[0][1][-1])))+ str(" RESIDUES HAD TO BE TRIMMED IN TOTAL TO CREATE THE FINAL PROTEIN. ********** \n" ))
            print("THESE ARE THE WINDOWS USED TO BUILD YOUR IDEALIZED COILED COIL PROTEIN:")
            while i< len(my[0][1]) - 1:
                window=str(my[0][1][i])
                print(window[10:-4])
                windows_for_suggested_protein.write(">"+str(window[10:-4])+str("\n"))
                i+=1
            windows_for_suggested_protein.close()
            list_of_windows_used=my[0][1][0:-1]
            number_of_proteins_made+=1
            initial_letter=1
            dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values={}
            




            #print(items[0])


























        #res = [key for key in dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values if dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values[key] == CC_prot_with_fewest_residues_trimmed] #LIST OF CC PROTS WITH FEWEST RESIDUES TRIMMED
        #print("THIS IS CC PROT WITH FEWEST RESIDUES TRIMMED")
        #CC_prot_with_fewest_residues_trimmed=res[0]
        #print(CC_prot_with_fewest_residues_trimmed)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
######IF REMAINDER IS 0 ###########################










if remainder ==0:
    print("REMAINDER IS 0.")
    
    last_window_used=initial_letter
    number_of_residues_trimmed=0
    iterator=0
    dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values= collections.OrderedDict()
#######ADD WHILE LOOP HERE TO BUILD N NUMBER OF PROTS. DROP FROM STORED_WINDOWS.
    number_of_proteins_made=0
    while number_of_proteins_made < number_of_proteins_to_build:
        
        
    
        for key in list_of_windows_used: #REMOVES THE WINDOWS WHICH WERE USED TO MAKE THE PREVIOUS PROTEIN
           
            for second_key in stored_a_windows:
                if key==second_key[0]:
                    stored_a_windows.remove(second_key)
            for second_key in stored_b_windows:
                if second_key[0]== key:
                    stored_b_windows.remove(second_key)
            for second_key in stored_c_windows:
                if second_key[0]== key:
                    stored_c_windows.remove(second_key)
            for second_key in stored_d_windows:
                if second_key[0]== key:
                    stored_d_windows.remove(second_key)
            for second_key in stored_e_windows:
                if second_key[0]== key:
                    stored_e_windows.remove(second_key)
            for second_key in stored_f_windows:
                if second_key[0]== key:
                    stored_f_windows.remove(second_key)
            for second_key in stored_g_windows:
                if second_key[0]== key:
                    stored_g_windows.remove(second_key)
               
        print("THESE ARE THE WINDOWS REMAINING AFTER THE PREVIOUS ITERATION")
        print("REMAINING A WINDOWS")
        print(stored_a_windows)
        print("REMAINING B WINDOWS")
        print(stored_b_windows)
        print("REMAINING C WINDOWS")
        print(stored_c_windows)
        print("REMAINING D WINDOWS")
        print(stored_d_windows)
        print("REMAINING E WINDOWS")
        print(stored_e_windows)
        print("REMAINING F WINDOWS")
        print(stored_f_windows)
        print("REMAINING G WINDOWS")
        print(stored_g_windows)
        #### MAYBE HERE
        if len(stored_a_windows)==0 and len(stored_b_windows)==0 and len(stored_c_windows)==0 and len(stored_d_windows)==0 and len(stored_e_windows)==0 and len(stored_f_windows)==0 and len(stored_g_windows)==0:
            print("THERE ARE NO REMAINING WINDOWS. ALL QUERY ANTIGENS WERE COVERED VIA THE LAST PROTEIN. EXITING NOW")
            sys.exit()
        else:
            pass
                
          
            
        while initial_letter <= 7: #USE THIS WHILE LOOP TO BUILD ALL 7 VARIANTS OF YOUR NTH IDEAL CC PROT 
            last_window_used = initial_letter
            a_windows=copy.deepcopy(stored_a_windows) #WINDOWS GET POPPED WITH EACH ITERATION TO BUILD THE POSSIBLE IDEAL CC. YOU NEED TO BE ABLE TO RESTORE THE LIST BY CALLING ON THE STORED LIST
            b_windows=copy.deepcopy(stored_b_windows)
            c_windows=copy.deepcopy(stored_c_windows)
            d_windows=copy.deepcopy(stored_d_windows)
            e_windows=copy.deepcopy(stored_e_windows)
            f_windows=copy.deepcopy(stored_f_windows)
            g_windows=copy.deepcopy(stored_g_windows)
            print("YEEHAW")
            print(iterator)
            print(a_windows)
            print(b_windows)
            print(c_windows)
            print(d_windows)
            print(e_windows)
            print(f_windows)
            print(g_windows)


        

            while iterator < number_of_windows_for_a_protein: ####MAKE SYS ARG
                ######MAYBE HERE 
                print("HERE YA GO")
                if len(stored_a_windows)==0 and len(stored_b_windows)==0 and len(stored_c_windows)==0 and len(stored_d_windows)==0 and len(stored_e_windows)==0 and len(stored_f_windows)==0 and len(stored_g_windows)==0:
                    
                #input("STOPPPP")
                if last_window_used== 1 or last_window_used== str(1) : #THE VALUE FOR A START WITH A AND END WITH A 
                    try: ### See if A list windows are empty WE ARE TRYING TO BUILD A PROTEIN WITH A AS THE FIRST REGISTER. TRY A LIST WINDOWS YOU TRIM NO RESIDUES
                        ideal_CC_prot_version_A.append(a_windows[0][1][0])
                        last_window_used=heptad_list[0]
                        fasta_header=str(a_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        a_windows.pop(0)
                        iterator+=1
                    except IndexError: ###  YOU NEED TO START WITH AN "A" LETTER.  "A" list windows are empty, try G list windows. You only have to trim 1 residue 
                        try:
                            current_frag=g_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[6]
                            fasta_header=str(g_windows[0][0])
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            g_windows.pop(0)
                            iterator+=1  
                        except IndexError: ### YOU NEED TO START WITH AN "A" LETTER. G list windows are empty, try F list windows. You trim 2 residues
                            try:
                                current_frag=f_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[5]
                                fasta_header=str(f_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                f_windows.pop(0)
                                iterator+=1
                            except IndexError: ### YOU NEED TO START WITH AN "A" LETTER. G AND F list windows are empty, try E list windows. You trim 3 residues
                                try:
                                    current_frag=e_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[4]
                                    fasta_header=str(e_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    e_windows.pop(0)
                                    iterator+=1
                                except IndexError: ###### G,F, AND, E list windows are empty, try D list windows. You trim 4 residues
                                    try:
                                        current_frag=d_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[3]
                                        fasta_header=str(d_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        d_windows.pop(0)
                                        iterator+=1
                                    except IndexError:
                                        try: ###### G,F,E,AND D list windows are empty, try C list windows You trim 5 residues
                                            current_frag=c_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag)
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[2]
                                            fasta_header=str(c_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            c_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ###### G,F,E,D,AND C list windows are empty, try B list windows You trim 6 residues
                                                current_frag=b_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[1]
                                                fasta_header=str(b_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                b_windows.pop(0)
                                                iterator+=1
                if last_window_used ==2 or last_window_used== str(2): ### WE JUST USED A "B" GET ANOTHER "B" 
                    try:
                        ideal_CC_prot_version_A.append(b_windows[0][1][0])
                        last_window_used=heptad_list[1]
                        fasta_header=str(b_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        b_windows.pop(0)
                        iterator+=1
                    except IndexError: ### B list windows are empty, try A list windows You trim 1 residue
                        try:
                            current_frag=a_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[0]
                            fasta_header=str(a_windows[0][0])
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            a_windows.pop(0)
                            iterator+=1
                        except IndexError: ### B AND A list windows are empty, try G list windows You trim 2 residues
                            try:
                                current_frag=g_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[6]
                                fasta_header=str(g_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                g_windows.pop(0)
                                iterator+=1
                            except IndexError: ###### B,A,AND G list windows are empty, try F list windows You trim 3 residues
                                try:
                                    current_frag=f_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[5]
                                    fasta_header=str(f_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    f_windows.pop(0)
                                    iterator+=1
                                except IndexError:
                                    try: ###### B,A,G,AND F list windows are empty, try E list windows. You trim 4 residues
                                        current_frag=e_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[4]
                                        fasta_header=str(e_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        e_windows.pop(0)
                                        iterator+=1
                                    except IndexError:
                                        try: ###### B,A,G,F AND E list windows are empty, try D list windows. You trim 5 residues
                                            current_frag=d_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag)
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[3]
                                            fasta_header=str(d_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            d_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ### B,A,G,F,E, AND D lists are empty, try C list You trim 6 residues
                                                current_frag=c_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[2]
                                                fasta_header=str(c_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                c_windows.pop(0)
                                                iterator+=1
                if last_window_used ==3 or last_window_used== str(3):
                    try:
                        ideal_CC_prot_version_A.append(c_windows[0][1][0])
                        last_window_used=heptad_list[2]
                        fasta_header=str(c_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        c_windows.pop(0)
                        iterator+=1
                    except IndexError: ### C list windows are empty, try B list windows You trim 1 residue
                        try:
                            current_frag=b_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[1]
                            fasta_header=str(b_windows[0][0])
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            b_windows.pop(0)
                            iterator+=1
                        except IndexError: ### C AND B list windows are empty, try A list windows You trim 2 residues
                            try:
                                current_frag=a_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[0]
                                fasta_header=str(a_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                a_windows.pop(0)
                                iterator+=1
                            except IndexError: ###### C,B,AND A list windows are empty, try G list windows You trim 3 residues
                                try:
                                    current_frag=g_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[6]
                                    fasta_header=str(g_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    g_windows.pop(0)
                                    iterator+=1
                                except IndexError:
                                    try: ###### C,B,A,AND G list windows are empty, try F list windows. You trim 4 residues
                                        current_frag=f_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[5]
                                        fasta_header=str(f_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        f_windows.pop(0)
                                        iterator+=1 
                                    except IndexError:
                                        try: ###### C,B,A,G, AND F list windows are empty, try E list windows. You trim 5 residues
                                            current_frag=e_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag)
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[4]
                                            fasta_header=str(e_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            e_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ### C,B,A,G,F, AND E lists are empty, try D list You trim 6 residues
                                                current_frag=d_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[3]
                                                fasta_header=str(d_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                d_windows.pop(0)
                                                iterator+=1


                if last_window_used ==4 or last_window_used== str(4): ### WE JUST USED A "D" GET ANOTHER "D"
                    try:
                        ideal_CC_prot_version_A.append(d_windows[0][1][0])
                        last_window_used=heptad_list[3]
                        fasta_header=str(d_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        d_windows.pop(0)
                        iterator+=1
                    except IndexError: ### D list windows are empty, try C list windows You trim 1 residue
                        try:
                            current_frag=c_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[2]
                            fasta_header=str(c_windows[0][0])
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            c_windows.pop(0)
                            iterator+=1
                        except IndexError: ### D AND C list windows are empty, try B list windows You trim 2 residues
                            try:
                                current_frag=b_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[1]
                                fasta_header=str(b_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                b_windows.pop(0)
                                iterator+=1
                            except IndexError: ###### D,C,AND B list windows are empty, try A list windows You trim 3 residues
                                try:
                                    current_frag=a_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[0]
                                    fasta_header=str(a_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    a_windows.pop(0)
                                    iterator+=1
                                except IndexError:
                                    try: ###### D,C,B,AND A list windows are empty, try G list windows. You trim 4 residues
                                        current_frag=g_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[6] 
                                        fasta_header=str(g_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        g_windows.pop(0)
                                        iterator+=1
                                    except IndexError:
                                        try: ###### D,C,B,A, AND G list windows are empty, try F list windows. You trim 5 residues
                                            current_frag=f_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag) 
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[5]
                                            fasta_header=str(f_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            f_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ### D,C,B,A,G, AND F lists are empty, try E list You trim 6 residues
                                                current_frag=e_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[4]
                                                fasta_header=str(e_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                e_windows.pop(0)
                                                iterator+=1
                if last_window_used ==5 or last_window_used== str(5): ### WE JUST USED A "E" GET ANOTHER "E"
                    try:
                        ideal_CC_prot_version_A.append(e_windows[0][1][0])
                        last_window_used=heptad_list[4]
                        fasta_header=str(e_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        e_windows.pop(0)
                        iterator+=1
                    except IndexError: ### E list windows are empty, try D list windows You trim 1 residue
                        try:
                            current_frag=d_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[3]
                            fasta_header=str(d_windows[0][0]) 
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            d_windows.pop(0)
                            iterator+=1
                        except IndexError: ### E AND D list windows are empty, try C list windows You trim 2 residues
                            try:
                                current_frag=c_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[2]
                                fasta_header=str(c_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                c_windows.pop(0)
                                iterator+=1
                            except IndexError: ###### E,D,AND C list windows are empty, try B list windows You trim 3 residues
                                try:
                                    current_frag=b_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[1]
                                    fasta_header=str(b_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    b_windows.pop(0)
                                    iterator+=1
                                except IndexError:
                                    try: ###### E,D,C,AND B list windows are empty, try A list windows. You trim 4 residues
                                        current_frag=a_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[0] 
                                        fasta_header=str(a_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        a_windows.pop(0)
                                        iterator+=1
                                    except IndexError:
                                        try: ###### E,D,C,B, AND A list windows are empty, try G list windows. You trim 5 residues
                                            current_frag=g_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag)
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[6]
                                            fasta_header=str(g_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            g_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ### E,D,C,B,A, AND G lists are empty, try F list You trim 6 residues
                                                current_frag=f_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[5]
                                                fasta_header=str(f_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                f_windows.pop(0)
                                                iterator+=1
                if last_window_used ==6 or last_window_used== str(6): ### WE JUST USED A "F" GET ANOTHER "F"
                    try:
                        ideal_CC_prot_version_A.append(f_windows[0][1][0])
                        last_window_used=heptad_list[5]
                        fasta_header=str(f_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        f_windows.pop(0)
                        iterator+=1
                    except IndexError: ### F list windows are empty, try E list windows You trim 1 residue
                        try:
                            current_frag=e_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[4]
                            fasta_header=str(e_windows[0][0])
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            e_windows.pop(0)
                            iterator+=1
                        except IndexError: ### F AND E list windows are empty, try D list windows You trim 2 residues
                            try:
                                current_frag=d_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[3]
                                fasta_header=str(d_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                d_windows.pop(0)
                                iterator+=1
                            except IndexError: ###### F,E,AND D list windows are empty, try C list windows You trim 3 residues
                                try:
                                    current_frag=c_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[2]
                                    fasta_header=str(c_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    c_windows.pop(0)
                                    iterator+=1
                                except IndexError:
                                    try: ###### F,E,D,AND C list windows are empty, try B list windows. You trim 4 residues
                                        current_frag=b_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[1]
                                        fasta_header=str(b_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        b_windows.pop(0)
                                        iterator+=1
                                    except IndexError:
                                        try: ###### F,E,D,C, AND B list windows are empty, try A list windows. You trim 5 residues
                                            current_frag=a_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag)
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[0]
                                            fasta_header=str(a_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            a_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ### F,E,D,C,B, AND A lists are empty, try G list You trim 6 residues
                                                current_frag=g_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[6]
                                                fasta_header=str(g_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                g_windows.pop(0)
                                                iterator+=1


                if last_window_used == 7 or last_window_used== str(7): ### WE JUST USED A "G" GET ANOTHER "G"
                    try:
                        ideal_CC_prot_version_A.append(g_windows[0][1][0])
                        last_window_used=heptad_list[6]
                        fasta_header=str(g_windows[0][0])
                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                        new_fasta_header=new_fasta_header[1]
                        new_fasta_header=new_fasta_header.replace(".out","")
                        new_fasta_header=str(">")+new_fasta_header
                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                        g_windows.pop(0)
                        iterator+=1
                    except IndexError: ### G list windows are empty, try F list windows You trim 1 residue
                        try:
                            current_frag=f_windows[0][1][0]
                            trimmed_frag=current_frag[1:]
                            ideal_CC_prot_version_A.append(trimmed_frag)
                            number_of_residues_trimmed+=1
                            last_window_used=heptad_list[5]
                            fasta_header=str(f_windows[0][0])
                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                            new_fasta_header=new_fasta_header[1]
                            new_fasta_header=new_fasta_header.replace(".out","")
                            new_fasta_header=str(">")+new_fasta_header
                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                            f_windows.pop(0)
                            iterator+=1
                        except IndexError: ### G AND F list windows are empty, try E list windows You trim 2 residues
                            try:
                                current_frag=e_windows[0][1][0]
                                trimmed_frag=current_frag[2:]
                                ideal_CC_prot_version_A.append(trimmed_frag)
                                number_of_residues_trimmed+=2
                                last_window_used=heptad_list[4]
                                fasta_header=str(e_windows[0][0])
                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                new_fasta_header=new_fasta_header[1]
                                new_fasta_header=new_fasta_header.replace(".out","")
                                new_fasta_header=str(">")+new_fasta_header
                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                e_windows.pop(0)
                                iterator+=1
                            except IndexError: ###### G,F,AND E list windows are empty, try D list windows You trim 3 residues
                                try:
                                    current_frag=d_windows[0][1][0]
                                    trimmed_frag=current_frag[3:]
                                    ideal_CC_prot_version_A.append(trimmed_frag)
                                    number_of_residues_trimmed+=3
                                    last_window_used=heptad_list[3]
                                    fasta_header=str(d_windows[0][0])
                                    new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                    new_fasta_header=new_fasta_header[1]
                                    new_fasta_header=new_fasta_header.replace(".out","")
                                    new_fasta_header=str(">")+new_fasta_header
                                    windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                    d_windows.pop(0)
                                    iterator+=1
                                except IndexError:
                                    try: ###### G,F,E,AND D list windows are empty, try C list windows. You trim 4 residues
                                        current_frag=c_windows[0][1][0]
                                        trimmed_frag=current_frag[4:]
                                        ideal_CC_prot_version_A.append(trimmed_frag)
                                        number_of_residues_trimmed+=4
                                        last_window_used=heptad_list[2]
                                        fasta_header=str(c_windows[0][0])
                                        new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                        new_fasta_header=new_fasta_header[1]
                                        new_fasta_header=new_fasta_header.replace(".out","")
                                        new_fasta_header=str(">")+new_fasta_header
                                        windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                        c_windows.pop(0)
                                        iterator+=1
                                    except IndexError:
                                        try: ###### G,F,E,D AND C list windows are empty, try B list windows. You trim 5 residues
                                            current_frag=b_windows[0][1][0]
                                            trimmed_frag=current_frag[5:]
                                            ideal_CC_prot_version_A.append(trimmed_frag)
                                            number_of_residues_trimmed+=5
                                            last_window_used=heptad_list[1]
                                            fasta_header=str(b_windows[0][0])
                                            new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                            new_fasta_header=new_fasta_header[1]
                                            new_fasta_header=new_fasta_header.replace(".out","")
                                            new_fasta_header=str(">")+new_fasta_header
                                            windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                            b_windows.pop(0)
                                            iterator+=1
                                        except IndexError: ### G,F,E,D,C, AND B lists are empty, try A list You trim 6 residues
                                                current_frag=c_windows[0][1][0]
                                                trimmed_frag=current_frag[6:]
                                                ideal_CC_prot_version_A.append(trimmed_frag)
                                                number_of_residues_trimmed+=6
                                                last_window_used=heptad_list[0]
                                                fasta_header=str(a_windows[0][0])
                                                new_fasta_header=fasta_header.split(("high_id_w_",1)[0])
                                                new_fasta_header=new_fasta_header[1]
                                                new_fasta_header=new_fasta_header.replace(".out","")
                                                new_fasta_header=str(">")+new_fasta_header
                                                windows_used_for_ideal_CC_prot_version_A.append(fasta_header)
                                                a_windows.pop(0)
                                                iterator+=1















            #nYAYYY
            final_ideal_CC_prot="".join(ideal_CC_prot_version_A)
            windows_used_for_ideal_CC_prot_version_A.append(float(number_of_residues_trimmed))

            dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values[str(final_ideal_CC_prot)]=windows_used_for_ideal_CC_prot_version_A
            ideal_CC_prot_version_A=[] #RESET TO BUILD CC PROT STARTING WITH WINDOW B, THEN WINDOW C, THEN WINDOW D...
            number_of_residues_trimmed=0 #RESET THE TRIM RESIDUE COUNT
            iterator=0 #RESET TO 0 SO YOU CAN GET THE K HEPTAD WINDOWS, FOR YOUR NTH CC PROTEIN 
            windows_used_for_ideal_CC_prot_version_A=[] #RESET TO GET WINDOWS USED TO BUILD CC PROT STARTING WITH WINDOW B, THEN WINDOW C, THEN WINDOW D...

            initial_letter+=1
        print("THIS IS THE DICT OF CC PROTS CREATED")
        print(dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values)
        #CC_prot_with_fewest_residues_trimmed = min(dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values.values()) 

        min_value = min(dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values.values(), key=lambda x: x[-1])[-1]
        dict_of_CC_proteins_w_fewest_residues_trimmed = {k: v for k, v in dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values.items() if v[-1] == min_value}
        print("THESE ARE THE PROTEINS WITH FEWEST RESIDUES TRIMMED")
        print(dict_of_CC_proteins_w_fewest_residues_trimmed)
        
        if len(dict_of_CC_proteins_w_fewest_residues_trimmed.keys()) >1:
            
            print("THERE IS MORE THAN ONE CC PROTEIN WHICH HAS THE FEWEST RESIDUES TRIMMED. WRITING THESE PROTEINS TO A TEMPORARY TEXT FILE AND RUNNING MARCOIL ANALYSIS TO SELECT THE PROTEIN WITH THE HIGHEST COILED COIL PROBABILITY.")
            temp_cc_prot_iterator=0 #SET TO 0 SO YOU CAN USE THE NUMBER IN THE FASTA HEADER TO EXTRACT BY KEY INDEX
            print(os.getcwd())
            os.chdir("/Users/jayspencer96/MARCOIL") #####MAKE SYS ARG
            print(os.getcwd())
            temporary_text_file_for_MARCOIL_analysis=open('temp.txt', 'w')
            for key in dict_of_CC_proteins_w_fewest_residues_trimmed:
                print("THIS IS KEY")
                print(key)
                temporary_text_file_for_MARCOIL_analysis.write(">CC_protein_" +str(temp_cc_prot_iterator)+str("\n"))
                temporary_text_file_for_MARCOIL_analysis.write(str(key) +str("\n"))
                temp_cc_prot_iterator+=1
            temporary_text_file_for_MARCOIL_analysis.close()

            os.system("./marcoil +dlsS temp.txt")
            os.chdir('/Users/jayspencer96/MARCOIL/Outputs/')###NOTEE output_for_MARCOIL must be an alias in your bash_profile that navigates to MARCOIL output
            print(os.getcwd())
            domains_file=open('Domains', 'r')




            dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values={}
            while True:

                        for line in domains_file:
                            if ">" in line:
                                current_CC_protein=str(line.strip())
                                list_of_domain_lengths=[]

                            elif "NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 2.0 :" in line:
                                number_of_2_percent_domains= int(line[-2]) ###Use this to determine how many lines you have to iterate through
                                if number_of_2_percent_domains == 0: ### No 2 percent domains. 

                                    list_of_domain_lengths.append(0)

                                elif number_of_2_percent_domains >=1:

                                    lines_to_iterate_through=1
                                    total_2_percent_domain_length=0
                                    while lines_to_iterate_through <= number_of_2_percent_domains:
                                        #print("THIS is LINE")
                                        #print(line)
                                        line=next(domains_file)
                                        #print("THIS IS NEXT LINE")
                                        #print(line)
                                        current_2_percent_domain_length=str(re.search(domain_length_pattern, line).group(0)) #GET THE NAME OF THE HIGH ID MATCH
                                        current_2_percent_domain_length = int(re.search(d_number_pattern, current_2_percent_domain_length).group(1))
                                        total_2_percent_domain_length=total_2_percent_domain_length+current_2_percent_domain_length
                                        lines_to_iterate_through+=1
                                    list_of_domain_lengths.append(total_2_percent_domain_length)



                            elif "NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 10.0 :" in line:
                                number_of_10_percent_domains= int(line[-2]) ###Use this to determine how many lines you have to iterate through

                                if number_of_10_percent_domains == 0: ### No 10 percent domains. 
                                    list_of_domain_lengths.append(0)


                                elif number_of_10_percent_domains >=1:
                                    lines_to_iterate_through=1
                                    total_10_percent_domain_length=0
                                    while lines_to_iterate_through <= number_of_10_percent_domains:
                                        line=next(domains_file)
                                        current_10_percent_domain_length=str(re.search(domain_length_pattern, line).group(0)) #GET THE NAME OF THE HIGH ID MATCH
                                        current_10_percent_domain_length = int(re.search(d_number_pattern, current_10_percent_domain_length).group(1))
                                        #current_10_percent_domain_length=int(current_10_percent_domain_length.group(0))
                                        total_10_percent_domain_length=total_10_percent_domain_length+current_10_percent_domain_length
                                        lines_to_iterate_through+=1
                                    list_of_domain_lengths.append(total_10_percent_domain_length)
                            elif "NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 50.0 :" in line:
                                number_of_50_percent_domains= int(line[-2]) ###Use this to determine how many lines you have to iterate through

                                if number_of_50_percent_domains == 0: ### No 50 percent domains. Go to the next line

                                    list_of_domain_lengths.append(0)


                                elif number_of_50_percent_domains >=1:
                                    lines_to_iterate_through=1
                                    total_50_percent_domain_length=0
                                    while lines_to_iterate_through <= number_of_50_percent_domains:
                                        line=next(domains_file)
                                        current_50_percent_domain_length=str(re.search(domain_length_pattern, line).group(0))  #GET THE NAME OF THE HIGH ID MATCH
                                        current_50_percent_domain_length = int(re.search(d_number_pattern, current_50_percent_domain_length).group(1))
                                        total_50_percent_domain_length=total_50_percent_domain_length+current_50_percent_domain_length
                                        lines_to_iterate_through+=1
                                    list_of_domain_lengths.append(total_50_percent_domain_length)
                            elif "NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 90.0 :" in line:
                                number_of_90_percent_domains= int(line[-2]) ###Use this to determine how many lines you have to iterate through

                                if number_of_90_percent_domains == 0: ### No 90 percent domains. Go to the next line
                                    list_of_domain_lengths.append(0)


                                elif number_of_90_percent_domains >=1:
                                    lines_to_iterate_through=1
                                    total_90_percent_domain_length=0
                                    while lines_to_iterate_through <= number_of_90_percent_domains:
                                        line=next(domains_file)
                                        current_90_percent_domain_length=str(re.search(domain_length_pattern, line).group(0))  #GET THE NAME OF THE HIGH ID MATCH
                                        current_90_percent_domain_length = int(re.search(d_number_pattern, current_90_percent_domain_length).group(1))
                                        total_90_percent_domain_length=total_90_percent_domain_length+current_90_percent_domain_length
                                        lines_to_iterate_through+=1
                                    list_of_domain_lengths.append(total_90_percent_domain_length)

                            elif "NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 99.0 :" in line:
                                number_of_99_percent_domains= int(line[-2]) ###Use this to determine how many lines you have to iterate through

                                if number_of_99_percent_domains == 0: ### No 10 percent domains. Go to the next line
                                    list_of_domain_lengths.append(0)
                                    dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values[str(current_CC_protein)] =list_of_domain_lengths

                                elif number_of_99_percent_domains >=1:
                                    lines_to_iterate_through=1
                                    total_99_percent_domain_length=0
                                    while lines_to_iterate_through <= number_of_99_percent_domains:
                                        line=next(domains_file)
                                        current_99_percent_domain_length=str(re.search(domain_length_pattern, line).group(0))  #GET THE NAME OF THE HIGH ID MATCH
                                        current_99_percent_domain_length = int(re.search(d_number_pattern, current_99_percent_domain_length).group(1))
                                        total_99_percent_domain_length=total_99_percent_domain_length+current_99_percent_domain_length
                                        lines_to_iterate_through+=1
                                    list_of_domain_lengths.append(total_99_percent_domain_length)
                                    dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values[str(current_CC_protein)] =list_of_domain_lengths

                            else:
                                pass
                        print("THIS IS VACCINE PROTS AND DOMAINS")
                        print(dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values)
                        break


            avg_dict = {k: (sum(v) / len(v)) for k, v in dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values.items()}

            print(avg_dict)
            best_protein=max(avg_dict, key=avg_dict.get)
            key_for_best_protein=int(int(best_protein[-1])-1)

            print("THIS IS AVG DICT")
            print(avg_dict)

            print("THIS PROTEIN HAS THE HIGHEST CC AVG SO IT WILL BE USED") #THERE MIGHT BE MORE THAN ONE, BUT AT THIS POINT YOU'VE COMPARED ENOUGH
            print(max(avg_dict, key=avg_dict.get))
            print("THIS IS THE INDEX NUMBER FOR THE BEST PROTEIN")
            print(key_for_best_protein)
            #items = list(ordered.items())
            #items = list(dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values.items())
            items = list(dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values.items())
            #FINAL_BEST_DICT=(items[key_for_best_protein])
            #print(items[dict_of_vaacine_proteins_as_keys_and_CC_domain_lengths_as_values[key_for_best_protein]])
            final_dict={} ### START FRESH SO THERE IS ONLY EVER 1 KEY IN THIS DICT, AND THAT'S THE BEST PROTEIN
            final_dict=items[key_for_best_protein]  #THIS WORKSSSSSS

            print("THIS IS THE BEST IDEALIZED COILED COIL PROTEIN")
            print(final_dict[0])
            suggested_proteins_file=open(PATH_TO_FINAL_OUTPUT+str("Suggested_proteins.txt"), "a+")
            print(len(final_dict[1]))
            i=0

            suggested_proteins_file.write(">Suggested_Coiled_Coil_Protein_"+str(count_for_suggested_proteins)+str("\n"))
            suggested_proteins_file.write(str(final_dict[0])+str("\n"))
            count_for_suggested_proteins+=1
            suggested_proteins_file.close()
            windows_for_suggested_protein=open(PATH_TO_FINAL_OUTPUT+str("Windows_used_to_make_suggested_proteins.txt"),"a+")
            windows_for_suggested_protein.write("********** THESE ARE THE WINDOWS USED TO MAKE SUGGESTED PROTEIN "+str(count_for_suggested_proteins) + str(". ")+str(int((final_dict[1][-1])))+ str(" RESIDUES HAD TO BE TRIMMED IN TOTAL TO CREATE THE FINAL PROTEIN. ********** \n" ))
            print("THESE ARE THE WINDOWS USED TO BUILD YOUR IDEALIZED COILED COIL PROTEIN:")
            while i< len(final_dict[1]) - 1:
                window=str(final_dict[1][i])
                print(window[10:-4])
                high_identity_file= open(str(MAIN_DIRECTORY) + str("HIGH_PERCENTAGE_IDENTITY/")+str(window), "r")
                for line in high_identity_file:
                    antigens_similar_to_vaccine.write(str(window[10:-4]) +str("\t") +str(line) +str("\n"))
                windows_for_suggested_protein.write(">"+str(window[10:-4])+str("\n"))
                i+=1
            windows_for_suggested_protein.close()
            list_of_windows_used=final_dict[1][:-1]
            #print("THIS IS LIST OF WINDOWS USED")
            print(list_of_windows_used)
            dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values={} #try clearing out here
            ideal_CC_prot_version_A=[]
            windows_used_for_ideal_CC_prot_version_A=[]
            number_of_residues_trimmed=0 #RESET THE TRIM RESIDUE COUNT
            iterator=0
            initial_letter=1
            last_window_used=1

            number_of_proteins_made+=1
        
        
        
        
        
        
        elif len(dict_of_CC_proteins_w_fewest_residues_trimmed.keys()) == 1:
            print("THERE IS ONLY ONE CC PROTEIN WHICH HAS THE FEWEST RESIDUES TRIMMED.")
            final_dict={} ### START FRESH SO THERE IS ONLY EVER 1 KEY IN THIS DICT, AND THAT'S THE BEST PROTEIN
            #final_dict=items[key_for_best_protein]  #THIS WORKSSSSSS

            print("THIS IS THE BEST IDEALIZED COILED COIL PROTEIN")
            my=list(dict_of_CC_proteins_w_fewest_residues_trimmed.items())
            print(my[0][0]) #THE LIST OF WINDOWS USED AND THE NUM RES TRIMMED AT END
            #print(my[0][1][i])
            
            suggested_proteins_file=open(PATH_TO_FINAL_OUTPUT+str("Suggested_proteins.txt"), "a+")
            
            i=0

            suggested_proteins_file.write(">Suggested_Coiled_Coil_Protein_"+str(count_for_suggested_proteins)+str("\n"))
            suggested_proteins_file.write(str(my[0][0])+str("\n"))
            count_for_suggested_proteins+=1
            suggested_proteins_file.close()
            windows_for_suggested_protein=open(PATH_TO_FINAL_OUTPUT+str("Windows_used_to_make_suggested_proteins.txt"),"a+")
            windows_for_suggested_protein.write("********** THESE ARE THE WINDOWS USED TO MAKE SUGGESTED PROTEIN "+str(count_for_suggested_proteins) + str(". ")+str(int((my[0][1][-1])))+ str(" RESIDUES HAD TO BE TRIMMED IN TOTAL TO CREATE THE FINAL PROTEIN. ********** \n" ))
            print("THESE ARE THE WINDOWS USED TO BUILD YOUR IDEALIZED COILED COIL PROTEIN:")
            

            while i< len(my[0][1]) - 1:
                window=str(my[0][1][i])
                print(window)
                print(window[10:-4])
                windows_for_suggested_protein.write(">"+str(window[10:-4])+str("\n"))
                high_identity_file= open(str(MAIN_DIRECTORY) + str("HIGH_PERCENTAGE_IDENTITY/")+str(window), "r")
                for line in high_identity_file:
                    antigens_similar_to_vaccine.write(str(window[10:-4]) +str("\t") +str(line) +str("\n"))
                i+=1
            windows_for_suggested_protein.close()
            list_of_windows_used=my[0][1][0:-1]
            initial_letter=1
            iterator=0
            last_window_used=1
            dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values={}
            number_of_proteins_made+=1
            
            
antigens_similar_to_vaccine.close()



            #print(items[0])


























        #res = [key for key in dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values if dict_of_CC_prots_seqs_as_keys_and_number_of_trimmed_residues_as_values[key] == CC_prot_with_fewest_residues_trimmed] #LIST OF CC PROTS WITH FEWEST RESIDUES TRIMMED
        #print("THIS IS CC PROT WITH FEWEST RESIDUES TRIMMED")
        #CC_prot_with_fewest_residues_trimmed=res[0]
        #print(CC_prot_with_fewest_residues_trimmed)
