#!/usr/bin/env python
import pickle, argparse
import site_frequency_spectrum as SFS

parser = argparse.ArgumentParser(description="This script iterates through a slim output file and for each simulation extracts the site frequency spectrum for each different element type and stores this information as a pickle jar.")

parser.add_argument("-i","--input", 
		required = True,
		dest = "input", 
		type =str, 
		help = "The sfs dict that has been generated using the script, get_SFS_from_slim.py")
parser.add_argument("-o","--output", 
		required = False,
		dest = "output", 
		type =str, 
		help = "The name of the SFS summary output file")

args = parser.parse_args()

slim_dict = pickle.load(open(args.input,"rb"))
organs = slim_dict[slim_dict.keys()[0]].keys()

## This time round, sfs_dict is not quite the same as before
## SFS dict is now a list of the slim runs (eg. 123_TEMP_SLIM.txt)
sel_dict = {}
neu_dict = {}

for slim in slim_dict.keys():
	for key in slim_dict[slim].keys():
		if key not in sel_dict.keys():
			neu_dict[key] = slim_dict[slim][key][0]
			sel_dict[key] = slim_dict[slim][key][1]
		else:
			neu_temp = SFS.merge_SFS( neu_dict[key] , slim_dict[slim][key][0] )
			sel_temp = SFS.merge_SFS( sel_dict[key] , slim_dict[slim][key][1] )
			neu_dict[key] = neu_temp
			sel_dict[key] = sel_temp

if not args.output:
	output_name = args.input.replace(".sfs.pkl",".summary.sfs")
else:
	output_name = args.output

output = open(output_name,"w")

for org in sorted(organs):
	sel_string = " ".join(map(str,sel_dict[org]))
	neu_string = " ".join(map(str,neu_dict[org]))
	output.write(org+"\n"+str(len(sel_dict[org]))+"\nselected sites\n"+sel_string+"\nneutral_sites\n"+neu_string+"\n")


