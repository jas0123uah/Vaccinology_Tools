""" This program parses MARCOIL results looking for "regular heptads". The heptad fragments must be at least a user-specified length or greater. The residues must also have a probabililty greater than or equal to a user-specified threshold.
1. ADD ABILITY TO LOOK AT RES NUMBER
Created by: Jay Spencer
 """
import re, sys, copy, os, shutil

####################################  PLACE THINGS USER MUST EDIT HERE ####################################  
in_fi=open(sys.argv[1], "r")
output=open(sys.argv[2]+str("/MARCOIL_CC_FRAGS.txt"), "w") 
No_CC=open(sys.argv[2]+str("/No_CC_test.txt"), "w")
heptad_output=open(sys.argv[2]+str("/MARCOIL_CC_FRAGS_HEPTAD_REGISTER.txt"), "w")
minimum_frag_size=sys.argv[3]
PATH_TO_MARCOIL=str("/Volumes/Mac_air/Mac_Air/MARCOIL/")
####################################  PLACE THINGS USER MUST EDIT ABOVE ####################################
#no_hash=find = re.compile(r"^(.*   ##)\..*")
percentage_pattern= re.compile(r"\d+\.\d")
sep = ' ' #eliminate the number from the ProbList FASTA HEADER

heptad_letter_pattern=re.compile(r"[a-g]\n")
AA_pattern=re.compile(r"[A-Z]")
heptad_list=["a\n","b\n","c\n","d\n", "e\n", "f\n", "g\n"]
number_for_initial_letter_of_frag=0
frag_count_for_current_seq=1 # You might experience breaks in "regular heptads" in the new antigen. Start at 1 for each new antigen 
current_frag_amino_acid_letter_list=[]
current_frag_heptad_letter_list=[]
heptad_init=0
test_if_input_has_any_CC=0 ###Sometimes the sequence will have no significant CC probability anywhere. The user needs to know this. 
current_seq= None




for line in in_fi:
    if ">" in line:
        if test_if_input_has_any_CC == 0 and current_seq != None: ### Make sure we are even looking at a sequence before entering this. If the previous antigen, (current_seq), had no CC frags we need to know. We also need to know there was a previous antigen (current_seq) to begin with. So we say current_seq != None
            #print("THIS IS THE CURRENT SEQ") ###The sequence with no CC frags.
            #print(current_seq)
            
            
            No_CC.write(current_seq)
            No_CC.write("\n")
            current_seq=str(line.strip()) #The last antigen did not have any CC frags based on current user-parameters. You're now at a new antigen in the ProbList file. Get its name.
            current_seq = current_seq.split(sep, 1)[0]
            frag_count_for_current_seq=1
            heptad_letter=None
            current_frag_heptad_letter_list=[] ### Frag is not regular or too short reset and try again
            current_frag_amino_acid_letter_list=[]
            number_for_initial_letter_of_frag=0 #maybe
            test_if_input_has_any_CC=0 #DIDNT DO ANYTHING
            heptad_init=0


        
        elif len(current_frag_amino_acid_letter_list)>=int(sys.argv[3]): ###This statement is useful because if you're parsing the next sequence's results you need to write out the last fragment for the previous sequence and reset everything.
                
                header_to_write= str(str(current_seq)+("_regular_fragment_")+str(frag_count_for_current_seq)+str("\n"))
                pep_seq="".join(current_frag_amino_acid_letter_list) #THE AMINO ACIDS FOR THE CURRENT FRAG JOINED TO A SINGLE STRING
                output.write(header_to_write)
                output.write(pep_seq)
                output.write("\n")
                h_letter_no_newline=str(current_frag_heptad_letter_list[0])
                h_letter_no_newline=h_letter_no_newline.strip()
                if h_letter_no_newline == "a":
                    heptad_output.write(str("1")+str("\n"))
                elif h_letter_no_newline == "b":
                    heptad_output.write(str("2")+str("\n"))
                elif h_letter_no_newline == "c":
                    heptad_output.write(str("3")+str("\n"))
                elif h_letter_no_newline == "d":
                    heptad_output.write(str("4") +str("\n"))
                elif h_letter_no_newline == "e":
                    heptad_output.write(str("5")+str("\n"))
                elif h_letter_no_newline == "f":
                    heptad_output.write(str("6")+str("\n"))
                elif h_letter_no_newline == "g":
                    heptad_output.write(str("7")+str("\n"))
                
                frag_count_for_current_seq=1 #### STARTING FRESH
                current_frag_heptad_letter_list=[]
                current_frag_amino_acid_letter_list=[]
                number_for_initial_letter_of_frag=0
                heptad_init=0
                test_if_input_has_any_CC=0
                current_seq=str(line.strip())
                current_seq = current_seq.split(sep, 1)[0]
        else:  # For the initial sequence in your MARCOIL results
            current_seq=str(line.strip())
            current_seq = current_seq.split(sep, 1)[0]
            
            #print("THIS IS CURRENT SEQ")
            #print(current_seq)
            frag_count_for_current_seq=1
            heptad_letter=None
            current_frag_heptad_letter_list=[] ### Frag is not regular or too short reset and try again
            current_frag_amino_acid_letter_list=[]
            number_for_initial_letter_of_frag=0 #maybe
            heptad_init=0
            test_if_input_has_any_CC=0
            


      
        
        
    elif re.search(percentage_pattern, line)!=None:  #If there are percentages in the line
        percentage=re.search(percentage_pattern, line)
        
        percentage=str(percentage.group(0))
        #print(percentage)
        
        previous_heptad_letter=heptad_letter ### Store the heptad letter from the previous iteration. From line 39, heptad_letter starts out equaling None.
        #print("THE PREVIOUS HEPTAD LETTER WAS")
        #print(previous_heptad_letter) ##Is NONE on the first iteration
        #print("THE CURRENT SEQ IS")
        #print(current_seq)
        #print("THE CURRENT LINE IS")
        #print(line)
        #print(current_frag_amino_acid_letter_list)
        if previous_heptad_letter!= None:
            index_for_initial_heptad_letter=heptad_list.index(previous_heptad_letter) ####Stores the index for the previous heptad letter so you can ensure a regular heptad repeat. .index returns the position at the first occurrence of the specified value
        else:
            index_for_initial_heptad_letter= None 
       
        heptad_letter=re.search(heptad_letter_pattern, line) #Gets the heptad letter for the current line
        
       
        heptad_letter=str(heptad_letter.group(0)) 
        
        #print("THIS IS HEPTAD LETTER")
        #print(heptad_letter)
        #print("THIS is index for previous heptad letter")
        #print(index_for_initial_heptad_letter)
        if index_for_initial_heptad_letter==6: ## Restarts the cycle once you get to 'g'
                index_for_initial_heptad_letter= -1

        if float(percentage) >= 2.0 and heptad_init ==0: ###We are trying to start a new heptad/fragment here
            current_frag_heptad_letter_list.append(heptad_letter)
            
            current_AA=re.search(AA_pattern, line)

            current_AA=str(current_AA.group(0))
           
            current_frag_amino_acid_letter_list.append(current_AA)  #LETTER 50 GETS APPENDED, BUT YOU NEVER WRITE
            index_for_initial_heptad_letter=heptad_list.index(heptad_letter)
            #print("THIS IS THE INDEX FOR THE INITIAL HEPTAD LETTER ")
            #print(index_for_initial_heptad_letter)
            if index_for_initial_heptad_letter==6: ## Restarts the cycle once you get to 'g'
                index_for_initial_heptad_letter= -1

            heptad_init+=1 ###Ensure this is only executed at the beginning of CC frags


        elif (float(percentage) >= 2.0 and heptad_letter == heptad_list[(int(index_for_initial_heptad_letter)+1)]): # You're in a regular heptad repeat
            current_AA=re.search(AA_pattern, line)
            #print("CHECKED THE LIST")

            current_AA=str(current_AA.group(0))
            #print("THIS IS THE CURRENT AA")
            #print(current_AA)
           
            current_frag_amino_acid_letter_list.append(current_AA)  
            
        else: # Either you are not in a regular repeat or the percent probability is below the threshold
            if len(current_frag_amino_acid_letter_list)< int(sys.argv[3]): #See how large the current fragment is. If it's small, clear everything out.
                #print("WE ARE HERRE")
                #print(current_seq)
                current_frag_heptad_letter_list=[] ### Frag is not regular or too short reset and try again
                current_frag_amino_acid_letter_list=[]
                heptad_init=0
                if float(percentage) >= 2.0: ### The current fragment was too small, but we need to see if the probability of the current line is at or above threshold. If it is we need to append the heptad letter and amino acid, as this might be the start of a regular heptad.
                    current_frag_heptad_letter_list.append(heptad_letter)
                    current_AA=re.search(AA_pattern, line)
                    current_AA=str(current_AA.group(0))
                    current_frag_amino_acid_letter_list.append(current_AA)
            else: #If the current fragment is big enough, write it to your output file.
            
                print(current_frag_amino_acid_letter_list)
                header_to_write= str(str(current_seq)+("_regular_fragment_")+str(frag_count_for_current_seq)+str("\n"))
                pep_seq="".join(current_frag_amino_acid_letter_list)
                output.write(header_to_write)
                output.write(pep_seq)
                output.write("\n")
                h_letter_no_newline=str(current_frag_heptad_letter_list[0])
                h_letter_no_newline=h_letter_no_newline.strip()
                if h_letter_no_newline == "a":
                    heptad_output.write(str("1")+str("\n"))
                elif h_letter_no_newline == "b":
                    heptad_output.write(str("2")+str("\n"))
                elif h_letter_no_newline == "c":
                    heptad_output.write(str("3")+str("\n"))
                elif h_letter_no_newline == "d":
                    heptad_output.write(str("4") +str("\n"))
                elif h_letter_no_newline == "e":
                    heptad_output.write(str("5")+str("\n"))
                elif h_letter_no_newline == "f":
                    heptad_output.write(str("6")+str("\n"))
                elif h_letter_no_newline == "g":
                    heptad_output.write(str("7")+str("\n"))
                test_if_input_has_any_CC=1 #FIXED?
                frag_count_for_current_seq+=1 #### incase you do find more frags in the current seq
                number_for_initial_letter_of_frag=0 #### Gets you back to the first 
                current_frag_heptad_letter_list=[]
                current_frag_amino_acid_letter_list=[]
                heptad_init=0
                
    else: ####GETS THE FINAL FRAGMENT IN THE MARCOIL RESULTS IF IT IS BIG ENOUGH
        if len(current_frag_amino_acid_letter_list)>=int(sys.argv[3]):
            
            header_to_write= str(str(current_seq)+("_regular_fragment_")+str(frag_count_for_current_seq)+str("\n"))
            pep_seq="".join(current_frag_amino_acid_letter_list)
            output.write(header_to_write)
            output.write(pep_seq)
            output.write("\n")
            h_letter_no_newline=str(current_frag_heptad_letter_list[0])
            h_letter_no_newline=h_letter_no_newline.strip()
            if h_letter_no_newline == "a":
                heptad_output.write(str("1")+str("\n"))
            elif h_letter_no_newline == "b":
                heptad_output.write(str("2")+str("\n"))
            elif h_letter_no_newline == "c":
                heptad_output.write(str("3")+str("\n"))
            elif h_letter_no_newline == "d":
                heptad_output.write(str("4") +str("\n"))
            elif h_letter_no_newline == "e":
                heptad_output.write(str("5")+str("\n"))
            elif h_letter_no_newline == "f":
                heptad_output.write(str("6")+str("\n"))
            elif h_letter_no_newline == "g":
                heptad_output.write(str("7")+str("\n"))
            
            #heptad_output.write(str(current_frag_heptad_letter_list[0]))
            #heptad_output.write("\n")
            frag_count_for_current_seq=1 #### STARTING FRESH
            test_if_input_has_any_CC=1 ######FIXED??
            current_frag_heptad_letter_list=[]
            current_frag_amino_acid_letter_list=[]
            number_for_initial_letter_of_frag=0
            heptad_init=0


   
output.close()
No_CC.close()
heptad_output.close()
equal_sized_frags=open(sys.argv[2]+str("/EQUAL_SIZED_FRAGS.txt"), "w") 
output=open(sys.argv[2]+str("/MARCOIL_CC_FRAGS.txt"), "r") 
for line in output:
    if ">" in line:
        header=line.strip()
        sequence=next(output)
        sequence=sequence.strip()
        if len(sequence) == int(sys.argv[3]):
            equal_sized_frags.write(header +str("\n") +str(sequence)+str("\n"))
        else:
            size_of_seq=len(sequence)
            iterator_for_size_of_seq=copy.deepcopy(int(sys.argv[3]))
            window_count=1
            start_point=0
            while iterator_for_size_of_seq <= size_of_seq:
                new_header=header+str("_window_")+str(window_count)
                new_seq=sequence[start_point:iterator_for_size_of_seq]
                equal_sized_frags.write(new_header +str("\n") +str(new_seq)+str("\n"))
                window_count+=1
                start_point+=1
                iterator_for_size_of_seq+=1
    else:
        pass
output.close()
equal_sized_frags.close()
path_to_equal_sized_frags=str(sys.argv[2]+str("/EQUAL_SIZED_FRAGS.txt"))

#####################DO MARCOIL ON EQUAL SIZED FRAGS TO ENSURE THEY ARE ALL REGULAR#########
os.chdir(PATH_TO_MARCOIL)
os.system(("./marcoil +dlsS " +str(path_to_equal_sized_frags)))
os.makedirs(str(sys.argv[2])+str("EQUAL_SIZED_FRAGS_MARCOIL_ANALYSIS/"), exist_ok=True)
shutil.move(str(PATH_TO_MARCOIL)+str("Outputs/Domains"), str(sys.argv[2])+str("EQUAL_SIZED_FRAGS_MARCOIL_ANALYSIS/Domains")) ## Used shutil bc Marcoil is on my external HD
shutil.move(str(PATH_TO_MARCOIL)+str("Outputs/ProbList"), str(sys.argv[2])+str("EQUAL_SIZED_FRAGS_MARCOIL_ANALYSIS/ProbList"))
shutil.move(str(PATH_TO_MARCOIL)+str("Outputs/ProbPerState"), str(sys.argv[2])+str("EQUAL_SIZED_FRAGS_MARCOIL_ANALYSIS/ProbPerState"))
           
                

################## PARSE MARCOIL RESULTS ON EQUAL SIZED FRAGS ###################                
            
irregular_windows=open(str(sys.argv[2])+str("EQUAL_SIZED_FRAGS_MARCOIL_ANALYSIS/IRREGULAR_HEPTAD_WINDOWS_AFTER_MAKING_WINDOWS_EQUAL_SIZE.txt"), "w")
number_for_initial_letter_of_frag=0
frag_count_for_current_seq=1 # You might experience breaks in "regular heptads" in the new antigen. Start at 1 for each new antigen 
current_frag_amino_acid_letter_list=[]
current_frag_heptad_letter_list=[]
heptad_init=0
test_if_input_has_any_CC=0 ###Sometimes the sequence will have no significant CC probability anywhere. The user needs to know this. 
current_seq= None  
equal_frag_file=open(str(sys.argv[2])+str("EQUAL_SIZED_FRAGS_MARCOIL_ANALYSIS/ProbList"), "r")
regular_equal_sized_frags=open(str(sys.argv[2])+str("EQUAL_SIZED_FRAGS_MARCOIL_ANALYSIS/REGULAR_WINDOWS.txt"), "w")
regular_equal_sized_frags_heptad_register=open(str(sys.argv[2])+str("EQUAL_SIZED_FRAGS_MARCOIL_ANALYSIS/REGISTER_FOR_REGULAR_WINDOWS.txt"), "w")


for line in equal_frag_file:
    if ">" in line:
        if test_if_input_has_any_CC == 0 and current_seq != None: ### Make sure we are even looking at a sequence before entering this. If the previous antigen, (current_seq), had no CC frags we need to know. We also need to know there was a previous antigen (current_seq) to begin with. So we say current_seq != None
            print("THIS IS THE CURRENT SEQ") ###The sequence with no CC frags.
            print(current_seq)
            
            print(line)
            
            
            irregular_windows.write(current_seq)
            irregular_windows.write("\n")
            current_seq=str(line.strip()) #The last fragment was irregular. You're now at a new fragment in the ProbList file. Get its name.
            current_seq = current_seq.split(sep, 1)[0]
            frag_count_for_current_seq=1
            heptad_letter=None
            current_frag_heptad_letter_list=[] ### Frag is not regular or too short reset and try again
            current_frag_amino_acid_letter_list=[]
            number_for_initial_letter_of_frag=0 #maybe
            test_if_input_has_any_CC=0 
            heptad_init=0


        
        elif len(current_frag_amino_acid_letter_list)==int(sys.argv[3]): ###This statement is useful because if you're parsing the next sequence's results you need to write out the last fragment for the previous sequence and reset everything. 
            header_to_write= str(str(current_seq)+("_regular_fragment_")+str(frag_count_for_current_seq)+str("\n"))
            pep_seq="".join(current_frag_amino_acid_letter_list) #THE AMINO ACIDS FOR THE CURRENT FRAG JOINED TO A SINGLE STRING
            regular_equal_sized_frags.write(header_to_write)
            regular_equal_sized_frags.write(pep_seq)
            regular_equal_sized_frags.write("\n")
            h_letter_no_newline=str(current_frag_heptad_letter_list[0])
            h_letter_no_newline=h_letter_no_newline.strip()
            if h_letter_no_newline == "a":
                regular_equal_sized_frags_heptad_register.write(str("1")+str("\n"))
            elif h_letter_no_newline == "b":
                heptad_output.write(str("2")+str("\n"))
            elif h_letter_no_newline == "c":
                regular_equal_sized_frags_heptad_register.write(str("3")+str("\n"))
            elif h_letter_no_newline == "d":
                regular_equal_sized_frags_heptad_register.write(str("4") +str("\n"))
            elif h_letter_no_newline == "e":
                regular_equal_sized_frags_heptad_register.write(str("5")+str("\n"))
            elif h_letter_no_newline == "f":
                regular_equal_sized_frags_heptad_register.write(str("6")+str("\n"))
            elif h_letter_no_newline == "g":
                regular_equal_sized_frags_heptad_register.write(str("7")+str("\n"))
            
            frag_count_for_current_seq=1 #### STARTING FRESH
            current_frag_heptad_letter_list=[]
            current_frag_amino_acid_letter_list=[]
            number_for_initial_letter_of_frag=0
            heptad_init=0
            
            #print(line)
            test_if_input_has_any_CC=0
            current_seq=str(line.strip())
            current_seq = current_seq.split(sep, 1)[0]
        else:  # For the initial sequence in your MARCOIL results
            current_seq=str(line.strip())
            current_seq = current_seq.split(sep, 1)[0]
           
            
            
            frag_count_for_current_seq=1
            heptad_letter=None
            current_frag_heptad_letter_list=[] ### Frag is not regular or too short reset and try again
            current_frag_amino_acid_letter_list=[]
            number_for_initial_letter_of_frag=0 #maybe
            heptad_init=0
            test_if_input_has_any_CC=0
            


      
        
        
    elif re.search(percentage_pattern, line)!=None:  #If there are percentages in the line
       
        percentage=re.search(percentage_pattern, line)
        
        percentage=str(percentage.group(0))
        print(percentage)
        
        previous_heptad_letter=heptad_letter ### Store the heptad letter from the previous iteration. From line 39, heptad_letter starts out equaling None.
        #print("THE PREVIOUS HEPTAD LETTER WAS")
        #print(previous_heptad_letter) ##Is NONE on the first iteration
        #print("THE CURRENT SEQ IS")
        #print(current_seq)
        #print("THE CURRENT LINE IS")
        #print(line)
        #print(current_frag_amino_acid_letter_list)
        if previous_heptad_letter!= None:
            index_for_initial_heptad_letter=heptad_list.index(previous_heptad_letter) ####Stores the index for the previous heptad letter so you can ensure a regular heptad repeat. .index returns the position at the first occurrence of the specified value
        else:
            index_for_initial_heptad_letter= None 
       
        heptad_letter=re.search(heptad_letter_pattern, line) #Gets the heptad letter for the current line
        
       
        heptad_letter=str(heptad_letter.group(0)) 
        
        #print("THIS IS HEPTAD LETTER")
        #print(heptad_letter)
        #print("THIS is index for previous heptad letter")
        #print(index_for_initial_heptad_letter)
        if index_for_initial_heptad_letter==6: ## Restarts the cycle once you get to 'g'
                index_for_initial_heptad_letter= -1

        if heptad_init ==0: ###We are trying to start a new heptad/fragment here
            current_frag_heptad_letter_list.append(heptad_letter)
            
            current_AA=re.search(AA_pattern, line)

            current_AA=str(current_AA.group(0))
           
            current_frag_amino_acid_letter_list.append(current_AA)  #LETTER 50 GETS APPENDED, BUT YOU NEVER WRITE
            index_for_initial_heptad_letter=heptad_list.index(heptad_letter)
            #print("THIS IS THE INDEX FOR THE INITIAL HEPTAD LETTER ")
            #print(index_for_initial_heptad_letter)
            if index_for_initial_heptad_letter==6: ## Restarts the cycle once you get to 'g'
                index_for_initial_heptad_letter= -1

            heptad_init+=1 ###Ensure this is only executed at the beginning of CC frags


        elif heptad_letter == heptad_list[(int(index_for_initial_heptad_letter)+1)]: # You're in a regular heptad repeat
            current_AA=re.search(AA_pattern, line)
            #print("CHECKED THE LIST")

            current_AA=str(current_AA.group(0))
            #print("THIS IS THE CURRENT AA")
            #print(current_AA)
           
            current_frag_amino_acid_letter_list.append(current_AA)  
            
        else: # Either you are not in a regular repeat or the percent probability is below the threshold
            if len(current_frag_amino_acid_letter_list)!= int(sys.argv[3]): #See how large the current fragment is. If it's small, clear everything out.
                #print("WE ARE HERRE")
                print(current_seq)
                current_frag_heptad_letter_list=[] ### Frag is not regular or too short reset and try again
                current_frag_amino_acid_letter_list=[]
                heptad_init=0
                """if float(percentage) >= 2.0: ### The current fragment was too small, but we need to see if the probability of the current line is at or above threshold. If it is we need to append the heptad letter and amino acid, as this might be the start of a regular heptad.
                    current_frag_heptad_letter_list.append(heptad_letter)
                    current_AA=re.search(AA_pattern, line)
                    current_AA=str(current_AA.group(0))
                    current_frag_amino_acid_letter_list.append(current_AA)"""
            else: #If the current fragment is the correct size write it to your output file.
            
                print(current_frag_amino_acid_letter_list)
                header_to_write= str(str(current_seq)+str("\n"))
                pep_seq="".join(current_frag_amino_acid_letter_list)
                regular_equal_sized_frags.write(header_to_write)
                regular_equal_sized_frags.write(pep_seq)
                regular_equal_sized_frags.write("\n")
                h_letter_no_newline=str(current_frag_heptad_letter_list[0])
                h_letter_no_newline=h_letter_no_newline.strip()
                if h_letter_no_newline == "a":
                    regular_equal_sized_frags_heptad_register.write(str("1")+str("\n"))
                elif h_letter_no_newline == "b":
                    regular_equal_sized_frags_heptad_register.write(str("2")+str("\n"))
                elif h_letter_no_newline == "c":
                    regular_equal_sized_frags_heptad_register.write(str("3")+str("\n"))
                elif h_letter_no_newline == "d":
                    regular_equal_sized_frags_heptad_register.write(str("4") +str("\n"))
                elif h_letter_no_newline == "e":
                    regular_equal_sized_frags_heptad_register.write(str("5")+str("\n"))
                elif h_letter_no_newline == "f":
                    regular_equal_sized_frags_heptad_register.write(str("6")+str("\n"))
                elif h_letter_no_newline == "g":
                    regular_equal_sized_frags_heptad_register.write(str("7")+str("\n"))
                test_if_input_has_any_CC=1 #FIXED?
                frag_count_for_current_seq+=1 #### incase you do find more frags in the current seq
                number_for_initial_letter_of_frag=0 #### Gets you back to the first 
                current_frag_heptad_letter_list=[]
                current_frag_amino_acid_letter_list=[]
                heptad_init=0
                
    else: ####GETS THE FINAL FRAGMENT IN THE MARCOIL RESULTS IF IT IS BIG ENOUGH
        if len(current_frag_amino_acid_letter_list)==int(sys.argv[3]):
            
            header_to_write= str(str(current_seq)+str("\n"))
            pep_seq="".join(current_frag_amino_acid_letter_list)
            regular_equal_sized_frags.write(header_to_write)
            regular_equal_sized_frags.write(pep_seq)
            regular_equal_sized_frags.write("\n")
            h_letter_no_newline=str(current_frag_heptad_letter_list[0])
            h_letter_no_newline=h_letter_no_newline.strip()
            if h_letter_no_newline == "a":
                regular_equal_sized_frags_heptad_register.write(str("1") +str("\n"))
            elif h_letter_no_newline == "b":
                regular_equal_sized_frags_heptad_register.write(str("2") +str("\n"))
            elif h_letter_no_newline == "c":
                regular_equal_sized_frags_heptad_register.write(str("3") +str("\n"))
            elif h_letter_no_newline == "d":
                regular_equal_sized_frags_heptad_register.write(str("4") +str("\n"))
            elif h_letter_no_newline == "e":
                regular_equal_sized_frags_heptad_register.write(str("5") +str("\n"))
            elif h_letter_no_newline == "f":
                regular_equal_sized_frags_heptad_register.write(str("6") +str("\n"))
            elif h_letter_no_newline == "g":
                regular_equal_sized_frags_heptad_register.write(str("7") +str("\n"))

            
            #heptad_output.write(str(current_frag_heptad_letter_list[0]))
            #heptad_output.write("\n")
            frag_count_for_current_seq=1 #### STARTING FRESH
            test_if_input_has_any_CC=1 ######FIXED??
            current_frag_heptad_letter_list=[]
            current_frag_amino_acid_letter_list=[]
            number_for_initial_letter_of_frag=0
            heptad_init=0       
       
equal_frag_file.close()
regular_equal_sized_frags.close()
regular_equal_sized_frags_heptad_register.close() 