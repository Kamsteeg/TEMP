#!/usr/bin/python
#Filename: CGM_Compare_GFF_V101.py
#Part of the compare gene model tool
#Maintainer Jeroen Kamsteeg

# python lib imports
import itertools
import time

class compare_GML:
	def __init__(self, gff3_main_dict, fasta_main_dict, ch_list, matchlvl):	
		################ temp vars
		self.search_range_type = 's'
		
		
		# imports
		self.gff3_main_dict = gff3_main_dict
		self.fasta_main_dict = fasta_main_dict
		self.ch_list = ch_list		
		
		# local vars
		self.gm_match_list = []	
		
		# start class functions		
		# -------------------------------------------		
		
		# compare the gene models 
		self.compareGML(matchlvl)
		
		# -------------------------------------------

	# compare the gff3 gene model items 
	def compareGML(self, matchlvl):	
		
		used_ref = []
		pos_ref = []
		
		# go tru the gene model dicts
		for ref, pre in itertools.combinations(self.gff3_main_dict, 2):
		#for ref, pre in itertools.permutations(self.gff3_main_dict, 2):
			
			used_ref.extend(str(ref))
			
			
			
			
					
			# loop true the referance gene model dict
			for ref_GM in range(1, len(self.gff3_main_dict[ref])):	
				if self.search_range_type == "ch": 
					# range is CH	
					if ref_GM in self.ch_list[ref]:				
						loc = self.ch_list[ref].index(ref_GM)				
						try:
							s_range = [self.ch_list[pre][loc], self.ch_list[pre][loc+1]]	
						except:
							pass
				if self.search_range_type == "s":
					if self.search_range_type == "s":
						range_nr = 50
					if self.search_range_type == "m":
						range_nr = 100
					if self.search_range_type == "l":
						range_nr = 250
					if ref_GM < range_nr: 
						s_range_start = 0
					if ref_GM > range_nr:
						s_range_start = ref_GM - range_nr
					if ref_GM + range_nr > len(self.gff3_main_dict[pre]):
						s_range_end = len(self.gff3_main_dict[pre])
					if ref_GM + range_nr < len(self.gff3_main_dict[pre]):
						s_range_end = ref_GM + range_nr
					s_range = [s_range_start, s_range_end]	
						
						
				# the size of the referance gene model
				ref_start = int(self.gff3_main_dict[ref][ref_GM][0][3])				
				ref_end = int(self.gff3_main_dict[ref][ref_GM][0][4])
				ref_size = ref_end - ref_start				
				# version info
				ref_version = self.gff3_main_dict[ref][ref_GM][0][0]			
				# get the number of exons
				ref_exon = self.gff3_main_dict[ref][ref_GM][1].count("exon") 
				# get the number of CDS 
				ref_cds = self.gff3_main_dict[ref][ref_GM][1].count("CDS")		
				# get the fasta id
				ref_gm_id = self.gff3_main_dict[ref][ref_GM][0][7]	
				
				#get fasta seq 				
				try:	
					ref_gm_seq = self.fasta_main_dict[ref][ref_gm_id].replace('*','')
				except:
					pass
					
				# var
				hit = 0
			
				# loop true the pre ch
				for pre_GM in range(s_range[0], s_range[1]):					
					try:					
						# size of the prediction gene model
						pre_start = int(self.gff3_main_dict[pre][pre_GM][0][3])
						pre_end = int(self.gff3_main_dict[pre][pre_GM][0][4])
						pre_size = pre_end - pre_start						
						# version info
						pre_version = self.gff3_main_dict[pre][pre_GM][0][0]
						# number of items in the gene model
						pre_len = len(self.gff3_main_dict[pre][pre_GM][1])					
						# get the number of exons
						pre_exon = self.gff3_main_dict[pre][pre_GM][1].count("exon")
						# get the number of CDS 
						pre_cds = self.gff3_main_dict[pre][pre_GM][1].count("CDS")						
						# get the fasta id 
						pre_gm_id = self.gff3_main_dict[pre][pre_GM][0][7]
						
						#get fasta seq 
						try:	
							pre_gm_seq = self.fasta_main_dict[pre][pre_gm_id].replace('*','')							
						except:
							pass
							
						# match type 1
						if int(matchlvl) == 1:	
							if ref_gm_seq == pre_gm_seq:
								# Add to the Hit dict
								item1 = str("gml_"+str(ref)+"_"+str(ref_GM))
								item2 = str("gml_"+str(pre)+"_"+str(pre_GM))		
								match_list = [item1, item2]										
								self.gm_match_list.extend([match_list])																					
								# control if there is a hit
								hit = 1
								# break
								break							 								
						
						# match type 2
						if int(matchlvl) == 2:
							findFASTA = ref_gm_seq.find(pre_gm_seq)
							if findFASTA == -1:
								pass
							else:
								# Add to the Hit dict
								item1 = str("gml_"+str(ref)+"_"+str(ref_GM))
								item2 = str("gml_"+str(pre)+"_"+str(pre_GM))		
								match_list = [item1, item2]										
								self.gm_match_list.extend([match_list])													
								# control if there is a hit
								hit = 1
								# break
								break							
						
						# match lvl 3
						if int(matchlvl) == 3:
							if pre_size == ref_size:
								findFASTA = ref_gm_seq.find(pre_gm_seq)
								if findFASTA == -1:
									pass
								else:
									# Add to the Hit dict
									item1 = str("gml_"+str(ref)+"_"+str(ref_GM))
									item2 = str("gml_"+str(pre)+"_"+str(pre_GM))		
									match_list = [item1, item2]										
									self.gm_match_list.extend([match_list])													
									# control if there is a hit
									hit = 1
									break							
						# match lvl 4
						if int(matchlvl) == 4:
							if pre_size == ref_size:
								if pre_exon == ref_exon	and pre_cds == ref_cds:
									item1 = str("gml_"+str(ref)+"_"+str(ref_GM))
									item2 = str("gml_"+str(pre)+"_"+str(pre_GM))		
									match_list = [item1, item2]										
									self.gm_match_list.extend([match_list])	
								break
						# match lvl 5
						if int(matchlvl) == 5:
							if pre_size == ref_size:
								if pre_len == ref_len:
									item1 = str("gml_"+str(ref)+"_"+str(ref_GM))
									item2 = str("gml_"+str(pre)+"_"+str(pre_GM))		
									match_list = [item1, item2]										
									self.gm_match_list.extend([match_list])	
								break  
						
													
					except:
						pass
						
						
				# add no match
				if hit == 0:
					item1 = str("gml_"+str(ref)+"_"+str(ref_GM))
					match_list = [item1]
					self.gm_match_list.extend([match_list])	
					
				# control if there is a hit
				hit = 0	
				
			
		for gml in self.gff3_main_dict:
			pos_ref.extend(gml)
			
		print pos_ref
		print used_ref
			
				
				
				
				
