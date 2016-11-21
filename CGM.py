#!/usr/bin/env python
#Filename: CGM.py
#Part of the compare gene model tool
#Maintainer Jeroen Kamsteeg

'''

Tool to combine gene model data from the same species  

'''

# import scripts/classes
from scripts.Load_Data import loadData
from scripts.Compare_GFF import compare_GML
from scripts.Combine_Results import combine_Results

# import functions, libs
import time
import sys
import argparse
import textwrap

def main(gff3_input, fasta_input, mlvl, outfile):
		
	# start class functions		
	# -----------------------------------------------------------
	
	print ""
	print "-----------------------------------------------------"
	print ""
			
	# step 1, load data 
	# -----------------------------------------------------------		
	
	startTime = time.time()
	print "STEP-1 Load data start"
	CGM_LD = loadData(gff3_input, fasta_input)		
	elapsedTime1 = time.time() - startTime
	print "STEP-1 Load finished ", int(elapsedTime1), "seconds"
	print ""

	# step 2, match gene models 
	# -----------------------------------------------------------
	
	
	
	
	startTime = time.time()
	print "STEP-2 Match Gene Models start"
	CGM_C = compare_GML(CGM_LD.gff3_main_dict, CGM_LD.fasta_main_dict, CGM_LD.ch_list, mlvl)
	elapsedTime2 = time.time() - startTime
	print "STEP-2 Match Gene Models finished ", int(elapsedTime2), "seconds"
	print ""

	# step 3, combine the match results
	# -----------------------------------------------------------
	
	startTime = time.time()
	print "STEP-3 Group the GMM start"
	CGM_CR = combine_Results(CGM_C.gm_match_list, CGM_LD.gff3_main_dict, CGM_LD.fasta_main_dict, outfile)
	elapsedTime3 = time.time() - startTime
	print "STEP-3 Group the GMM finished ", int(elapsedTime3), "seconds"
	
	print ""
	print "-----------------------------------------------------"

	# -----------------------------------------------------------



if __name__ == '__main__':	
	parser = argparse.ArgumentParser(prog='PROG', 
				formatter_class=argparse.RawDescriptionHelpFormatter,
				description=textwrap.dedent('''\
				Please do not mess up this text!
				--------------------------------
				I have indented it
				exactly the way
				I want it
				'''), epilog=textwrap.dedent('''
				--------------------------------					
				'''))
	parser.add_argument('--gff3', nargs='+', default=0, type = argparse.FileType('r'),
						help='add gff3 files to compare.                            ')
	parser.add_argument('--fasta', nargs='+', default=0, type = argparse.FileType('r'),
						help='add the coresponding fasta files.                     ')
	parser.add_argument('--match', default=1,
						help=textwrap.dedent('''
						the level of similaraties that are needed to match gene models : 1= size match, 2= size, cds, exon match 3= size pre sequense found in ref sequense, 4= size, presequence = refsequense match 
						'''))
	parser.add_argument('--out', default="cgm_results.txt",
						help=textwrap.dedent('''
						name of result file
						 '''))
	
	
	######
	
	arg = parser.parse_args()		
	main(arg.gff3, arg.fasta, arg.match, arg.out)
	
	"""
	try:	
		arg = parser.parse_args()		
		main(arg.gff3, arg.fasta, arg.match, arg.out)
	except:
		parser.print_help()
		sys.exit(0)

	"""
