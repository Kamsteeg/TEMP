

import time

class combine_Results:
	def __init__(self, GMML, GFFD, FASTAD, OUT):
		
		################ temp vars
		
		# name out file
		self.NOUTF = OUT
		
		# start class functions		
		# -------------------------------------------
		
		
		
		# group the gene model matches
		GMGL = self.group_results(GMML, GFFD)
		
		# make the output file
		self.OUTF = self.result_file_make()
		
		# orders the results 		 	
		self.order_results(GMGL, GFFD, FASTAD)
		
		# close the output file
		self.result_file_close()

		# -------------------------------------------
	
	# group the gene model matches
	# GMML = gene model match list
	def group_results(self, GMML, GFFD):
		
		
		
		
		GMGL = []
		
		for GMM in range(len(GMML)):
			NGMG = 1
			if GMML[GMM] in GMGL:
				NGMG = 0
			else:
				try:
					for GMG in range(len(GMGL)):
						if len(GMML[GMM]) == 1:
							if GMML[GMM] in GMGL[GMG]:
								NGMG = 0
								GMTL = GMML[GMM] + GMGL[GMG]			
								GMGL[GMG] = list(set(GMTL))
								break
						if len(GMML[GMM]) > 1:
							if GMML[GMM][0] in GMGL[GMG] or GMML[GMM][1] in GMGL[GMG]:
								NGMG = 0
								GMTL = GMML[GMM] + GMGL[GMG]
								GMGL[GMG] = list(set(GMTL))
								break					
				except:
					pass
					
			if NGMG == 1:
				GMGL.extend([GMML[GMM]])
				
		return GMGL
	
	
	# print the results
	# GMGL = gene model group list, GFFD = GFF dict, FASTAD = Fasta dict							 
	def order_results(self, GMGL, GFFD, FASTAD):
		
		
		
		for GMG in GMGL:
			FASTAL = []
			
			GMG_ID = ""
			
			#self.OUTF.write("GM_ID=")	
			#self.OUTF.write(str(GMG_ID))
			self.OUTF.write("\tGMD_ID=")			
			
			for GM in sorted(GMG):
				GM_RAW = GM.rstrip('\n').split('_')
				self.OUTF.write(str(GM))
				self.OUTF.write(" ")
				try:
					GM_info = GFFD[int(GM_RAW[1])][int(GM_RAW[2])]
					FASTAL.extend([FASTAD[int(GM_RAW[1])][GM_info[0][7]]])
				except:
					pass
			
			#FASTAL = sorted(FASTAL, key=len)
			
			#self.OUTF.write("\tREFFASTA=")	
			#self.OUTF.write(str(FASTAL[0]))					
			self.OUTF.write("\n")		
		
			
	# makes a file for the results
	def result_file_make(self):
		out_file = open(self.NOUTF, "w")
		return out_file
	
	# closes the result file	
	def result_file_close(self):
		self.OUTF.close()
				
		
	def result_print(self):
		pass

				
			



	
		
