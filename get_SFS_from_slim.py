#!/usr/bin/env python

import sys, argparse, gzip,tom_slim,pickle
from multiprocessing import Pool
from site_frequency_spectrum import SFS_from_frequencies as SFS
from tom import in_range, brace

def SFS_by_organ(raw_slim):
	x = tom_slim.slim(raw_slim)
	org_dict = x.organ_mutations()
	org_lengths = x.organ_lengths()
	sfs_dict = {}
	for j in org_dict.keys():
		freq_sel = []
		freq_neu = []
		for mut in org_dict[j]:
			if float(mut[3]) != 0.0:
				freq_sel.append(int(mut[7]))
			else:
				freq_neu.append(int(mut[7]))
		sel_sites = sel_prop_dict[j] * org_lengths[j]
		neu_sites = (1.0-sel_prop_dict[j]) * org_lengths[j]
		## Get the SFS for neutral and selected sites  using the lengths
		## of the different genomic elements
		sel_sfs = SFS(freq_sel,sel_sites,x.sampleN)
		neu_sfs = SFS(freq_neu,neu_sites,x.sampleN)

		sfs_dict[j] = [neu_sfs,sel_sfs]
	return [x.name, sfs_dict]	

parser = argparse.ArgumentParser(description="This script iterates through a slim output file and for each simulation extracts the site frequency spectrum for each different element type and stores this information as a pickle jar.")

parser.add_argument("-i","--input", 
		required = True,
		dest = "input", 
		type =str, 
		help = "The SLiM output file")
parser.add_argument("-o","--output", 
		required = False,
		dest = "output", 
		type =str, 
		help = "The name of the output pickle jar")
parser.add_argument("--gz", 
		required = False,
		dest = "gz", 
		action = "store_true", 
		help = "Is the input file gzipped? It is recommended (by me) that the file should be gzipped",
		default = False)
parser.add_argument("--procs","-p", 
		required = False,
		dest = "procs", 
		type = int, 
		help = "How many cores do you want to engage in this operation?",
		default = 1)

sel_prop_dict = {"g0":0.0,"g1":0.75,"g2":0.0,"g3":0.0} ## This is a dict of the proportion of mutations that are selected in a certain genomic element. Could get automate this from SLiM object...

args = parser.parse_args()


## Get the default output name if not supplied
if not args.output:
	if args.gz:
		output_name = args.input.replace(".txt.gz",".sfs.pkl")
	else:
		output_name = args.input.replace(".txt",".sfs.pkl")
else:
	output_name = args.output


## This is a new idea for me, if the user (me in the future, wants to extact the SFS for a gzipped set of slim runs, they can.  This little bed gets the right function for the job
if args.gz:
	slim_iterator = tom_slim.slim_reader_gzip
else:
	slim_iterator = tom_slim.slim_reader


if args.procs == 1:
	output_dict = {}
	for i in slim_iterator(args.input):
		output_list = SFS_by_organ(i)
		output_dict[output_list[0]] = output_list[1]
elif args.procs > 1:
	p = Pool(args.procs)
	output_dict = {}
	results = p.map(SFS_by_organ,slim_iterator(args.input))
	for i in results:
		output_dict[i[0]] = i[1]

pickle.dump( output_dict, open(output_name, "wb" ) )

