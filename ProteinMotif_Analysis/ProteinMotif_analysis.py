#!/usr/local/bin/python3
print("---------------------------------------------------------------------------------------------------------------------------")
print("\tThis Program can search, fetch and check the protein sequences from NCBI, then do some analysis as follows:")
print("\t\t(1) Calculate and Plot the Similarity of Protein Sequences\n\t\t(2) Blast Sequence Alignment\n\t\t(3) Generate a Phylogenetic Tree\n\t\t(4) Scan Protein Sequence of Interest with Motifs from the PROSITE Database")
print("---------------------------------------------------------------------------------------------------------------------------\n")
####Q1--------------------------------------------------------------------------------------------------------------------
import os, subprocess, re
#os.mkdir("protein_seq_analysis")
#os.system("cd protein_seq_analysis/")
###define the function to search and fetch sequences from NCBI

##define a function (precise-search) that need input specific protein name and taxonomy. 
def precise_get_seq(pr,taxo):
	input_echo="echo Your protein family is {} and your taxonomic subset is {}".format(pr,taxo)
	subprocess.call(input_echo, shell=True)	
	#Ask the user whether they want to get partial sequences or non-partial sequences?	
	if input("\nCHECK-1-Do you want to get all available sequences or non-partial sequences(recommended)?\nAll data,type: all\tNon-partial: just press 'Enter'\n")=="all":
		esearch_input=(pr+"[Protein Name] AND "+taxo+"[Organism]")
		print("Precise all-available search:")	
	else :
		print("Precise non-partial search:")
		esearch_input=(pr+"[Protein Name] AND "+taxo+"[Organism] NOT PARTIAL")
	esearch_command="esearch -db protein -query '{}'".format(esearch_input)
	subprocess.call(esearch_command, shell=True)
	return esearch_input

##define another function (open-search) that need input protein family and specific taxonomy.
def get_seq(pr,taxo):
	input_echo="echo Your protein family is {} and your taxonomic subset is {}".format(pr,taxo)
	subprocess.call(input_echo, shell=True)
	if input("\nCHECK-1-Do you want to get all available sequences or non-partial sequences(recommended)?\nAll data,type: all\tNon-partial: just press 'Enter'\n")=="all":
		esearch_input=(pr+" AND "+taxo+"[Organism]")
		print("Open all-available search:")
	else :
		esearch_input=(pr+" AND "+taxo+"[Organism] NOT PARTIAL")
		print("Open non-partial search:")
	esearch_command="esearch -db protein -query '{}'".format(esearch_input)
	subprocess.call(esearch_command, shell=True)
	
	return esearch_input

##define a function to fetch the sequences from NCBI
def fetch_seq(esearch_input):
	print("\nPlease wait, fetching the sequences you want...\n")
	efetch_command="esearch -db protein -query '{}' | efetch -format fasta > protein.fasta".format(esearch_input)
	subprocess.call(efetch_command, shell=True)



###while loop for sequence number check before fetching
##Use some starting values to let while loop get started.
while_start = True
while while_start:
##let the user to define the data set.
	pr=input("What is the protein family(default value is glucose-6-phosphatase)\n") or "glucose-6-phosphatase"
	taxo=input("What is the subset of the taxonomic tree(default value is Aves)\n") or "Aves"
##Choose precise search or open search?
##Get the sequences data and put into "protein.fasta"
	precise_or_not=input("\nPlease choose whether precise search or open search for protein do you want?\nprecise, tpye:p \t open, type:o\n")
	if precise_or_not.lower() == "p":
		print("Precise search:")
		esearch_input=f"{precise_get_seq(pr,taxo)}"
		print(esearch_input)
	elif precise_or_not.lower() == "o":
		print("Open search:")
		esearch_input=f"{get_seq(pr,taxo)}"
		print(esearch_input)
	else :
		print("I don't understand your choice, but just use the precise search.\n")
		esearch_input=f"{precise_get_seq(pr,taxo)}"
		print("Terms for esearch:\t",esearch_input)
##Check the sequence number
	print("\nCHECK-2-Sequence Number:")
	continue_reset = input("\nIn above search results,<Count> means sequence number you got. \nIf you got more than 1000 sequences, you will wait for a long time for sequence fetching.\nIf you want to continue,type:continue \nIf you want to reset the search information, type: reset\n")
	if continue_reset.lower()=="continue":
		while_start2 = True
	elif continue_reset.lower()=="reset":
		while_start2 = False
		print("Reset the search information...\n")
	else :
		print("I don't know,but maybe you want to continue.Let's continue.\n")
		while_start2= True
##While loop for species check
	while while_start2:
		fetch_seq(esearch_input)
		print("CHECK-3-Species Number:")
##open and read the fasta file and get meaningful column
		pr_fasta=open("protein.fasta")
		pr_seq=pr_fasta.read().split(">")
		pr_seq=pr_seq[1:]
##extract the species name into a list and set the list to calculate how many different species we get.
		species_name=[]
		for seq in pr_seq:
			match = re.search("\[(.+)\]",seq).group() #Extract [species name]
			species_name.append(match)   #Put [species name] into a list species_name
		species_set = list(set(species_name))  #Set the values in the list to get species kinds.
		print(f"There are {len(pr_seq)} sequences of {len(species_set)} species in your corrent dataset\n")
##Ask the user continue with the dataset or not?
		continue_or_not=input("Do you want to continue with the current dataset?\nType:continue \t Type: reset\n")
		if continue_or_not.lower()=="continue":
			while_start2 = False
			while_start = False
		elif continue_or_not.lower()=="reset":
			print("Reset the search information...\n")
			while_start = True
			break
		else :
			print("I don't know,but maybe you want to continue.Let's continue.\n")
			while_start2 = False
			while_start = False

##Print out the prepared results of data set.
print("Good job! Your starting dataset is prepared in file 'protein.fasta'.")


####Q2-----------------------------------------------------------------------------------------------------------------

###Choose the longest sequence of each species for further analysis
print("\nThis programm will automatically choose the longest sequence of each species in your starting data base for further analysis.\n")
##set a starting dictionary to save sequences using species as index.
species_seq_dic={}
for spec in species_set:
        species_seq_dic[spec] = ""

##put the longest sequence of each species into the dictionary.
for seq in pr_seq:
        match = re.search("\[(.+)\]",seq).group()
        seq = seq.replace(" ","_") #delate the accesstion number and use "_" to connect the protein and specie
        if len(seq) > len(species_seq_dic[match]):
                species_seq_dic[match] = seq
##put the sequences into a list
pr_seq_choose=list(species_seq_dic.values())

##get a check value for filtering.
import numpy
##calculat the median sequence length of all the sequences.And let the user choose which part to get.
pr_seq_choose_len = []
for seq in pr_seq_choose:
        pr_seq_choose_len.append(len(seq))
median = numpy.median(pr_seq_choose_len)
default = median - 30
print("The median number of sequence length is ", median)

check_value = int(input("You can input a length limit to remove those short sequences from your data set.\n Or use the defaults[median-30].\n") or default)
if check_value not in range(0,1000):
        check_value = defaults
print("Your check value is",check_value)

##extract the sequences in "pr_seq_choose_len" with good length.
print("\nThe length of sequences of choice are as follows.\n")
pr_sequence_choose = []
new_len = []
for seq in pr_seq_choose:
	if len(seq) >= check_value:
		new_len.append(len(seq))
		pr_sequence_choose.append(seq)
print(new_len)
output_str =">"+">".join(pr_sequence_choose)
choose_output = open("pr_seq_choose.fasta", "w")
choose_output.write(output_str)
print(f"There are {len(pr_sequence_choose)} output sequences in 'pr_seq_choose.fasta'")
choose_output.close()

print("\nWell done! Let's start analysis!\n")


####Q2.2-------------------------------------------------------------------------------------
##Level of Conservation
next_anal=input("\nDo you want to get the similarity and level of conservation of the protein sequences in the dataset?\nYes or No:\n") or "Yes"
if next_anal.lower()=="yes":
	print("Analysis-1-Calculate and Plot the Similarity of Protein Sequences:")
	plotcon_command1 = "plotcon pr_seq_choose.fasta -graph data"
	plotcon_command2 = "plotcon pr_seq_choose.fasta -graph x11"
	subprocess.call(plotcon_command1, shell=True)
	subprocess.call(plotcon_command2, shell=True)
	print("The conservation level data is in file 'plotcon1.dat'.")
else :
	print("OK, sikp this step.\n")

##Blast Sequence Alignment
next_anal=input("\nDo you want to do Blast sequence alignment within the dataset?\nYes or No:\n") or "Yes"
if next_anal.lower()=="yes":
	print("Analysis-2-Blast Sequence Alignment:")
	#create a database using our choose sequence in "pr_seq_choose.fasta"
	blast_db_command="makeblastdb -in pr_seq_choose.fasta -dbtype prot"
	subprocess.call(blast_db_command,shell=True)
	#blast the protein sequences using the database just generated
	blastp_command="blastp -query pr_seq_choose.fasta -db pr_seq_choose.fasta -out pr_seq_choose_blast.tsv -evalue 1e-5 -outfmt 7"
	subprocess.call(blastp_command ,shell=True)
	#show the blast result on the screen
	subprocess.call("less pr_seq_choose_blast.tsv", shell=True)
	print("The output file of Blast within your dataset is 'pr_seq_choose_blast.tsv', including:\nquery acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score")
else :
	print("OK, sikp this step.\n")

##phylogeny Tree
next_anal=input("\nDo you want to generate a Phylogenetic Tree using the current dataset?\nYes or No:\n") or "Yes"
if next_anal.lower()=="yes":
	print("Analysis-3-Generate a Phylogenetic Tree with Edialign:")
	edialign_command = "edialign -sequences pr_seq_choose.fasta -outfile pr_seq_choose.edialign -outseq edialign_pr_seq.fasta"
	subprocess.call(edialign_command, shell=True)
	subprocess.call("less pr_seq_choose.edialign",shell=True)
	print("The alignment sequence is in file 'edialign_pr_seq.fasta'")
	print("The Average sequence Length, Similarity, Phylogenetic Tree are in file 'pr_seq_choose.edialign'")
else :
	print("OK, sikp this step.\n")

##Scan protein sequence(s) of interest with motifs from the PROSITE database
next_anal=input("\nDo you want to scan a protein sequence of interest with motifs?\nYes or No:\n") or "Yes"
if next_anal.lower()=="yes":
	all_or_partial=input("You want to use all sequences in our dataset or choose some of the sequences of your interest(Default is partial)?\nType:all\tType:partial\n")
	if all_or_partial.lower() == "all":
		subprocess.call("rm -f pr_seq.patmatmotifs",shell=True) # clear the content of "pr_seq.patmatmotifs" in case of the error caused by executing this program for multiple times.
		print("Analysis-4-Scan Protein Sequence of Interest with Motifs from the PROSITE Database:")
		for seq in pr_sequence_choose:
			one_pr_seq = open("one_pr_seq.fasta","w")
			one_pr_seq.write(">" + seq)
			one_pr_seq.close()
		#scan the input protein sequence into outputfile"one_pr_seq.patmatmotifs", append the content of this file into "pr_seq.patmatmotifs".	
			motif_scan_command="patmatmotifs -sequence one_pr_seq.fasta -outfile one_pr_seq.patmatmotifs -full"
			subprocess.call(motif_scan_command,shell=True)
			subprocess.call("cat one_pr_seq.patmatmotifs >> pr_seq.patmatmotifs", shell=True)
	else :
		subprocess.call("rm -f pr_seq.patmatmotifs",shell=True) # clear the content of "pr_seq.patmatmotifs" in case of the error caused by executing this program for multiple times.
		print("Analysis-4-Scan Protein Sequence of Interest with Motifs from the PROSITE Database:")
		subprocess.call('grep ">" pr_seq_choose.fasta',shell=True) #Put the sequence information on the screen for reference.
		print("\nFrom all above sequence information, you can find the accession number.\n")
		while_again = True
		input_seq_list=[]
		while while_again == True:
			one_pr_seq = open("one_pr_seq.fasta","w")
			seq_info = input("Please choose one sequence that you're interested in and input the accession number.The default is XP_031456423.1\n") or "XP_031456423.1"
			for seq in pr_sequence_choose:
				if seq_info in seq:
					print("Find sequence:\n", seq)
					input_seq_list.append(seq)
					one_pr_seq.write(">" + seq)
			one_pr_seq.close()
			#scan the input protein sequence into outputfile"one_pr_seq.patmatmotifs", append the content of this file into "pr_seq.patmatmotifs".
			motif_scan_command="patmatmotifs -sequence one_pr_seq.fasta -outfile one_pr_seq.patmatmotifs -full"
			subprocess.call(motif_scan_command,shell=True)
			subprocess.call("cat one_pr_seq.patmatmotifs >> pr_seq.patmatmotifs", shell=True)
			again_or_not=input("\nDo you have other proteins to scan (Default is Yes)?\nYes or No?\n") or "Yes"
			if again_or_not.lower() == "yes":
				while_again = True
			else :
				while_again = False
				print("\nOk, protein sequence motif scanning finished.")
				print("Your input sequences are as follows:")
				print("\n".join(input_seq_list))
	print("Output results for all the protein you scanned are in file 'pr_seq.patmatmotifs'")
	subprocess.call("less pr_seq.patmatmotifs",shell=True)
	#Show the motif name:
	print("Name of motif(s) associated with the protein you interested in:")
	subprocess.call('grep "Motif = " pr_seq.patmatmotifs > motif_name.txt',shell=True)
	motif_name_txt=open("motif_name.txt")
	motif_name=motif_name_txt.read().split("\n")
	motif_name_list=list(set(motif_name[:-1]))
	print("All kinds of motifs found are in the following list:\n",motif_name_list)
	motif_name_txt.close()	
else :
        print("OK, sikp this step.\n")
pr_fasta.close()

