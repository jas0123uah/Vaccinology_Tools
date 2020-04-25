import sys, linecache
from Bio import SeqIO
perfect_cc_seqs=open(sys.argv[1], "r") #The input file you need to condense.
perfect_register=sys.argv[2] #The register corresponding to your input file.
path_to_output=str(sys.argv[3]) #The directory to write your output to.
counter=0
headers=[]
seqs_used=[]
registers=[]
with open(sys.argv[1], "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
            seq_register_list=[]
            counter+=1
            header=str(record.id)
            sequence=str(record.seq)
            if sequence in seqs_used:
                index_of_interest=seqs_used.index(sequence)
                original_header=str(headers[index_of_interest])
                new_header=original_header+str("_IDENTICAL_")+str(header)
                headers=[new_header if x==original_header else x for x in headers]
            elif sequence not in seqs_used:
                brand_new_header=str(">")+str(header)
                headers.append(brand_new_header)
                seqs_used.append(sequence)
                register=str(linecache.getline(sys.argv[2], counter))
                register.strip()
                registers.append(register)
condensed_seqs=dict(zip(headers, seqs_used))
condensed_seq_file=open(path_to_output+str("CONDENSED_PERFECT_CC_SEQS"), "w+")
condensed_register_file=open(path_to_output+str("REGISTER_FOR_CONDENSED_PERFECT_CC_SEQS"), "w+")
for key, value in condensed_seqs.items():
        condensed_seq_file.write(str(key)+str("\n")+str(value)+str("\n"))
condensed_seq_file.close()
for register in registers:
        condensed_register_file.write(register)
condensed_register_file.close()