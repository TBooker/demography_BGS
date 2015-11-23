import sys, argparse, gzip,tom_slim,tom,pickle
from multiprocessing import Pool
from site_frequency_spectrum import SFS_from_all_frequencies as SFS
def get_mutations(input):
	x = tom_slim.slim(i)
	print x.mutations

def in_range(point,range):
	if point >= range[0] and point <= range[1]:
		return True
	else:
		return False

#p = Pool(3)
#p.map(get_mutations, tom_slim.slim_reader_gzip(sys.argv[1]))
output_dict = {}
for i in tom_slim.slim_reader_gzip(sys.argv[1]):
	x = tom_slim.slim(i)
	org_dict = {}
	org_lengths = {}
	for mut in x.mutations:
	#	print mut
		for org in x.organs:
			if org[0] not in org_lengths.keys():
				org_lengths[org[0]]=int(org[2]-,int(org[1])+1
			else:
				org_lengths[org[0]]+=int(org[2]-,int(org[1])+1
			
			if in_range(int(mut[2]),[int(org[1]),int(org[2])]):
				if org[0] not in org_dict.keys():
					org_dict[org[0]] = [mut]
				else:
					org_dict[org[0]].append(mut)
	print org_dict.keys()	
	sfs_dict = {}
	
	for j in org_dict.keys():
		freq_sel = []
		freq_neu = []
		for mut in org_dict[j]:
			print mut
			if float(mut[3]) != 0.0:
				print mut
				freq_sel.append(int(mut[7]))
			else:
				freq_neu.append(int(mut[7]))
		sel_sfs = SFS(freq_sel,x.sampleN)
		neu_sfs = SFS(freq_neu,x.sampleN)
		sfs_dict[j] = [neu_sfs,sel_sfs]
	output_dict[x.name] = sfs_dict
	print output_dict
	break
