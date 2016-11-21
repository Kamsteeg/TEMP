#!/usr/bin/python
#Filename: CGM_Load_Data.py
#Part of the compare gene model tool
#Maintainer Jeroen Kamsteeg

#imports
import itertools

class loadData:
	
	def __init__(self, gff3_input_files, fasta_input_files):
		
		# start class functions		
		# -------------------------------------------						
	
		# makes all the needed variabels
		self.setVariables(gff3_input_files, fasta_input_files)

		# loop tru the input files and put them in a dictionary
		for file_id in self.gff3_main_dict:
			self.fillDictionaryGFF(gff3_input_files, file_id)
			self.fillDictionaryFasta(fasta_input_files, file_id)			
		
		# make the ch list
		for dic in range(len(self.gff3_main_dict)):			
			for item in range(1, (len(self.gff3_main_dict[dic]))):
				if item in self.ch_list[dic]:
					if self.search_range == len(self.ch_list[dic])-1:						
						range_start = self.ch_list[dic][self.search_range]
						range_stop = len(self.gff3_main_dict[dic]) 			
					else:			
						try:							
							range_start = self.ch_list[dic][self.search_range]
							range_stop = self.ch_list[dic][self.search_range+1]  			
							self.search_range = self.search_range + 1		
						except:
							pass
		
			self.search_range = 0
			
		# -------------------------------------------	

	# set al the variables and dicts, lists, etz
	def setVariables(self, gff3_input, fasta_input):
		#make dictionarys and lists based on the number of input files
		
		# GFF3 dict
		self.gff3_main_dict = { new_dict:{} for new_dict in range(len(gff3_input))}	
			
		# Fasta dict
		self.fasta_main_dict = { new_dict:{} for new_dict in range(len(gff3_input))}
			
		# ch range list
		self.ch_list = [[] for _ in range(len(gff3_input))]	
		
		self.search_range = 0
				
	# open a file
	def openFiles(self, infile):
		
		return infile
	
	# gets the data from the gff3 files	
	def fillDictionaryGFF(self, gff3_input_files, item_id):
		
		# open the gff file
		gff3file = gff3_input_files[item_id]			
		
		# loop trough gff3 file
		dict_key_count = 1
		ch_list_count = 1
		
		# gene model items
		new_gene_model = 0
				
		for line in gff3file:			
			if line[0] == "#":
				if line[0:17] == "##sequence-region":					
					ch_info = line.rstrip('\n').split(' ')
					self.ch_list[item_id].extend([dict_key_count])
				else:
					pass
			else:				
				# parse the gff3 data (gff3 is tab seperated)
				split_line = line.rstrip('\n').split('\t')				
								
				# fill the dictionary				
				if split_line[2].upper() == "GENE" and new_gene_model == 1:
					itemInfo = split_line[8].rstrip('\n').replace('=',';').replace(':',';').split(';')
					for item in range(len(itemInfo)):
						if itemInfo[item] == "ID":
							itemID = itemInfo[item + 2]	
							itemID = itemID.split('.')
							itemID = itemID[0]		
					
					#combine gene model items
					itemList = split_line[0:7]
					itemList.append(itemID)
					
					self.gff3_main_dict[item_id][dict_key_count] = [itemList, new_gene_model_item]
					dict_key_count = dict_key_count + 1
					new_gene_model = 0
					
				if split_line[2].upper() == "GENE" and new_gene_model == 0:
					new_gene_model = 1
					new_gene_model_item = []
				
				if split_line[2].upper() != "GENE" and new_gene_model == 1:	
					new_gene_model_item.append(split_line[2])
							
	# gets the data from the fasta files					
	def fillDictionaryFasta(self, fasta_input_files, item_id):
		fastaFile = fasta_input_files[item_id]
		for line in fastaFile:
			if line[0] == ">":
				split_line = line.rstrip('\n').replace('.',' ').split(' ')	
				header = split_line[0][1:len(split_line[0])]
				nextline = "seq"
			if line[0] != ">" and nextline == "seq":
				seq = line
				self.fasta_main_dict[item_id][header] = line.rstrip('\n')
				nextline = ""
				
				
				
				
				
				
