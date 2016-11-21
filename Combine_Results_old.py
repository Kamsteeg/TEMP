#!/usr/bin/python
#Filename: Combine_Results.py
#Part of the compare gene model tool
#Maintainer Jeroen Kamsteeg

class combine_Results:
	#
	def __init__(self, gmml, gffd, fastadict, outfile):
		
		################ temp vars
		
		# name out file
		self.outfile = outfile
		
		# start class functions		
		# -------------------------------------------
				
		# group the gene model matches
		gmgl = self.group_results(gmml)		
		# make the output file
		self.out = self.result_file_make()		
		# orders the results 		 	
		self.order_results(gmgl, gffd, fastadict)		
		# close the output file
		self.result_file_close()

		# -------------------------------------------
	
	# group the gene model matches
	# gmml = gene model match list
	def group_results(self, gmml):	
		gmgl = []
		for gmm in range(len(gmml)):
			ngmg = 1
			if gmml[gmm] in gmgl:
				ngmg = 0
			else:
				try:
					for gmg in range(len(gmgl)):
						if len(gmml[gmm]) == 1:
							if gmml[gmm] in gmgl[gmg]:
								ngmg = 0
								gmtl = gmml[gmm] + gmgl[gmg]			
								gmgl[gmg] = list(set(gmtl))
								break
						if len(gmml[gmm]) > 1:
							if gmml[gmm][0] in gmgl[gmg] or gmml[gmm][1] in gmgl[GMG]:
								ngmg = 0
								gmtl = gmml[gmm] + gmgl[gmg]
								gmgl[gmg] = list(set(gmtl))
								break					
				except:
					pass					
			if ngmg == 1:
				gmgl.extend([gmml[gmm]])				
		return gmgl
	
	
	# print the results
	# gmgl = gene model group list, gffd = GFF dict, fastadict = Fasta dict							 
	def order_results(self, gmgl, gffd, fastadict):		
		for gmg in gmgl:
			# new fasta list 
			FASTAL = []	
			# new start list
			startlist = []
			
			#####################################
			# temp name 
			gmg_id = ""
			####################################
						
			#self.out.write("GM_ID=")	
			#self.out.write(str(gmg_id))
			#self.out.write("\tGMD_ID=")						
			for GM in sorted(gmg):
				GM_RAW = GM.rstrip('\n').split('_')
				self.out.write(str(GM))
				self.out.write(" ")
				try:
					GM_info = gffd[int(GM_RAW[1])][int(GM_RAW[2])]
					
					FASTAL.extend([fastadict[int(GM_RAW[1])][GM_info[0][7]]])
					startlist.extend([GM_info[0][3]])
				except:
					pass
			#self.out.write("\tgm positions=")	
			#self.out.write(str(startlist))	
			#FASTAL = sorted(FASTAL, key=len)			
			#self.out.write("\tREFFASTA=")	
			#self.out.write(str(FASTAL[0]))					
			self.out.write("\n")				
			
	# makes a file for the results
	def result_file_make(self):
		out_file = open(self.outfile, "w")
		return out_file
	
	# closes the result file	
	def result_file_close(self):
		self.out.close()
