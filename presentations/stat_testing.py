#!/usr/bin/env python

# sel_testing is a testing suite to thoroughly test the function and accuracy of all of the
# python scripts used in my selection detection analysis pipeline
############################################################################################################################
from time import clock

from os import chdir as cd
from os import listdir as ls
from os import path
from cPickle import load
from time import clock

import sim_tools
import sim_parser

import argparse as ap
import bz2
import leslie
import numpy as np
import pop_classes
import pop_gen
import re
import subprocess as sp
import sys

import traceback

npar = np.array
divide = np.divide
np.seterr(divide = "ignore", invalid = "ignore") # don't show a warning for every zero division or infinite division error you get
np.set_printoptions(threshold = 100) # change the number of elements that triggers array summarization

if "phased" in sys.argv:
	PHASED = True
else:
	PHASED = False
if "unphased" in sys.argv:
	UNPHASED = True
else:
	UNPHASED = False
############################################################################################################################
start = clock()
############################################################################################################################
## Run a simple simulation with few variants, strong selection at center,
## high recombination, and two subpopulations
#param_file = "/Volumes/Data/2013-analysis-reboot/selection-detection/experiments/function-testing/function-testing-parameters.txt"
#with open(param_file, "r") as f:
#	params = f.read().strip()
#params = params.split()
## test the function of sim_launcher
#these_args = sim_parser.SimParser.parse_args(params)
#verbose = these_args.verbose
#del these_args.verbose
## check the simulation parameters
#sim_parser.SimCheck(these_args)
## print any necessary warnings
#sim_parser.SimWarnings(these_args)
## make a Sim object
#this_sim = sim_tools.Sim(these_args)
#print this_sim
## run the simulation
#print "***Starting simulations!", (clock()-start)/60.0, " min.***"
#this_sim.run()
#print "***Finished with simulations!", (clock()-start)/60.0, " min.***"
#print this_sim.filename
## -------------------------------------------
#print "-" * 20
#print "Simulations run, parsed, and filtered"
#print "sim_launcher.py ALL TESTS PASSED"
#print "-" * 20
############################################################################################################################
#DIR = this_sim.filename
DIR = "/Volumes/Data/2013-analysis-reboot/selection-detection/experiments/function-testing/test-suite_11"
#PREFIX = this_sim.filename
PREFIX = "test-suite_11"
cd(DIR)
# count the number of reps you find here
rep_re = re.compile(r"rep\.\d{3}\.pkl\.bz2") # make a pattern to match the rep files to
file_list = ls(".") # list all of the files in this directory
rep_files = [el for el in file_list if re.search(rep_re, el)] # keep the files that match the pattern
rep_count = len(rep_files) # count how many files you found
# unpickle the sim object
simfile = bz2.BZ2File(PREFIX+".Sim.pkl.bz2", "r")
this_sim = load(simfile)
# check to see if the number of repetitions matches between the files and the sim parameters
if rep_count != this_sim.n_repetitions:
	print """There are %d reps found, but this_sim says it needs %d reps;
	\ncontinuing with analysis of the %d found reps... ...""" % (rep_count, this_sim.n_repetitions, rep_count)

# scale the window size
args_window_size = [10000]
SCALED_WINDOW_SIZE = [10000.0 / 500000.0]

# Get the pickled sim and parsed popsets
repN = 0
repfile = bz2.BZ2File(PREFIX + ".rep" + ".%03d" % repN + ".pkl.bz2")
my_popset = load(repfile)
print my_popset
# make the popset a new popset, since I've probably changed the PopSet's functions
string_chromosomes = npar(my_popset.chromosomes, dtype="str")
# if selection is on, restore the selected site to the chromosomes array and make it into a string array
if my_popset.sel_position is not None:
	string_chromosomes[:, my_popset.sel_index] = my_popset.sel_genotypes

# now make a new popset object to use
my_popset = pop_classes.PhasedPopSet(chromosomes = string_chromosomes,
									 n_segregating_sites = my_popset.n_segregating_sites,
									 pop_indexes = my_popset.pop_indexes,
									 pop_names = my_popset.pop_names,
									 pop_sizes = my_popset.pop_sizes,
									 positions = my_popset.decimal_positions,
									 sel_position = my_popset.sel_position,
									 sel_pops = my_popset.sel_pops,
									 sel_origins = my_popset.sel_origins,
									 chromosome_indexes =my_popset.chromosome_indexes
									 )

my_popset.remove_invariant_sites()
print "Popset with invariant sites removed:"
print my_popset
my_popset.reinit(this_sim.seq_length) # daf_cutoff is default
# -------------------------------------------
print "-" * 20
print "One parsed rep successfully read back in"
print "parsed rep reading ALL TESTS PASSED"
print "-" * 20
############################################################################################################################
def tester(command):
	try:
		exec command in globals()
		print "\t", "-" * 20
		print "\tCommand Successful:", command
	except:
		print "\tCOMMAND FAILED:", command
		info = sys.exc_info()
		for file, lineno, function, text in traceback.extract_tb(info[2]):
			print "\t", file, "line", lineno, "in", function
			print "\t=>", repr(text)
		print "\t** %s: %s" % info[:2]

############################################################################################################################
# pick out some of the sites to look at
test_site_indexes = [my_popset.sel_index,
					 my_popset.common_site_indexes[round(my_popset.n_common_sites * 0.25)],
					 my_popset.common_site_indexes[round(my_popset.n_common_sites * 0.9)]]
test_site_indexes = npar(test_site_indexes)
############################################################################################################################
if PHASED:
	# open the haplogrouping log
	hgroup_log = open(PREFIX + ".haplogroup-log.txt", "w")
	print "PHASED POPSET TESTING"
	print my_popset.__str__().replace("\n", "\n\t")
	# Run each PhasedPopSet function
	#my_popset.make_ihhs()						#my_popset.ihh_ancestral[this_pop] = [list of ihhA per site]
	##											#my_popset.ihh_derived[this_pop] = [list of ihhD per site]
	tester("my_popset.make_hamming_matrix()")				#my_popset.hamming_matrix[site, pop1, pop2]
	# print my_popset.make_matrix_text(range(5, 20), range(25, 100))
	tester("my_popset.make_pis()")						#my_popset.pi_within[this_pop] = [list of pis by site]
	tester("my_popset.make_singletons()")					#my_popset.singletons[this_pop] = [list of singletons per site]
	tester("my_popset.make_derived_counts()")				#my_popset.derived_counts[this_pop] = [list of derived counts per site]
	#calculate things important only for more than 1 pop
	if len(my_popset.pop_names) > 1:
		print "\tcalculating for multiple pops", (clock() - start) / 60.0, " min.***"
		tester("my_deltadafs = my_popset.delta_dafs()")
		tester("my_global_fst = my_popset.global_fst()")
		tester("my_fsts = my_popset.fsts(this_sim.pop_pairs)")
		#tester("my_xpehhs = my_popset.xpehhs(this_sim.pop_pairs))
	##ihh based stats
	#tester("my_ihss = my_popset.iHSs()")
	#tester("my_delta_ihhs = my_popset.delta_ihhs()")
	tester("file_tuple = my_popset.make_rehh_input(this_sim.filename, rep_num=0)")
	tester("my_popset.run_rehh_script(file_tuple)")

	# these dictionaries will hold all of the values that are calculated for every site in the loops below
	all_dxys = {}				# window: [list of {pop: Dxy} at each site]
	all_d12s = {}				# window: [list of {pop: D12} at each site]
	all_pis = {}				# window: [list of {pop: pi} at each site]
	all_tajima_ds = {}			# window: [list of {pop: Tajima's D} at each site]
	all_fstars = {}				# window: [list of {pop: F*} at each site]
	all_dstars = {}				# window: [list of {pop: D*} at each site]
	all_hs = {}					# window: [list of {pop: (theta_H, fay_wu_H)} at each site]
	windowed_haplogroups = {}	# window: [list of Haplogroups objects for each site]
	for w in args_window_size:
		windowed_haplogroups[w] = []	#will be a list of Haplogroups objects, one for each site
		all_dxys[w] = []				#will be a list of Dxy dictionaries, one for each site, consisting of pop:siteDxy pairs
		all_d12s[w] = []				#will be a list of D12 dictionaries, one for each site, consisting of pop:siteD12 pairs

	print "\tcalculating haplogrouped stats", (clock() - start) / 60.0, " min.***"
	#classify all of the haplogroups and get all window indexes for later use
	for (posN, dpos, bpos) in zip(my_popset.common_site_indexes, my_popset.decimal_positions[my_popset.common_site_indexes], my_popset.bp_positions[my_popset.common_site_indexes]):
		# posN, dpos, bpos = my_popset.common_site_indexes[5], my_popset.decimal_positions[my_popset.common_site_indexes][5], my_popset.bp_positions[my_popset.common_site_indexes][5]
		print >> hgroup_log, "//\n" + "rep:%d\tpos_index:%d\tdecimal_pos:%f\tbp_pos:%d" % (repN, posN, dpos, bpos)
		# make a popset filtered for the derived allele at this site (use for haplogrouping etc.)
		tester("derived_filter = my_popset.get_people_with_allele_at_site(posN, pop_gen.DERIVED_ALLELE)")
		print >> hgroup_log, derived_filter
		#[/later] start with largest window and then cut it down from there, only calculating once?
		# loop over the window sizes you're using
		for (scaledW, unscaledW) in zip(SCALED_WINDOW_SIZE, args_window_size):
			# scaledW, unscaledW = SCALED_WINDOW_SIZE[0], args_window_size[0]
			#get the positions to keep for the window desired
			tester("kp, ki = my_popset.get_window_sites(dpos, scaledW)")
			#print >> hgroup_log, "kp: ", kp
			#print >> hgroup_log, "ki: ", ki
			#classify the haplogroups according to this site and this window
			tester("my_haplogroups = my_popset.get_haplogroups(ki, people_indexes=derived_filter, haplogroup_log=hgroup_log)")
			print >> hgroup_log, my_haplogroups
			# add this to the list of Haplogroups objects
			tester("windowed_haplogroups[unscaledW].append(my_haplogroups)")
			# calculate Dxy in the derived popset using the haplogrouped indices and window_indexes
			tester("all_dxys[unscaledW].append(my_popset.window_dxys(windowed_haplogroups[unscaledW][-1]))")
			# calculate D12 in the whole popset using the haplogrouped indices and window keepIs (in the hapset)
			tester("all_d12s[unscaledW].append(my_popset.window_d12s(windowed_haplogroups[unscaledW][-1]))")

	print "\tcalculating Garud et al. stats", (clock() - start) / 60.0, " min.***"
	my_garud_stats = {}
	for (scaledW, unscaledW) in zip(SCALED_WINDOW_SIZE, args_window_size):
		tester("my_garud_stats[unscaledW] = my_popset.garud_stats(scaledW)")

	print "\tcalculating windowed statistics", (clock() - start) / 60.0, " min.***"
	# haplogroup based stats
	for w in args_window_size:
		#these list comprehensions make a list of all the sites' stat values in one window (like doing a for loop over sites)
		#allele spectrum based stats
		tester("all_pis[w] = [my_popset.window_pi(hg.window_indexes) for hg in windowed_haplogroups[w]]")
		tester("all_tajima_ds[w] = [my_popset.window_tajima_d(hg.window_indexes) for hg in windowed_haplogroups[w]]")
		tester("all_dstars[w] = [my_popset.window_d_star(hg.window_indexes) for hg in windowed_haplogroups[w]]")
		tester("all_fstars[w] = [my_popset.window_f_star(hg.window_indexes) for hg in windowed_haplogroups[w]]")
		tester("all_hs[w] = [my_popset.window_fay_wu_h(hg.window_indexes) for hg in windowed_haplogroups[w]]")

	# test the printing function
	my_stat_dict = {"global.fst" : my_global_fst,
					"fst" : my_fsts,
					"delta.daf" : my_deltadafs,
					"hgroups" : windowed_haplogroups,
					"d12" : all_d12s,
					"d.star" : all_dstars,
					"dxy" : all_dxys,
					"h.output" : all_hs,
					"f.star" : all_fstars,
					"tajima.d" : all_tajima_ds,
					"window.pi" : all_pis,
					"garud" : my_garud_stats}
	tester("my_popset.save_stats(PREFIX, this_sim.pop_pairs, 0, my_stat_dict, args_window_size)")
############################################################################################################################
if UNPHASED:
	# open the haplogrouping log
	hgroup_log = open("unphased/" + PREFIX + ".haplogroup-log.txt", "w")
	print "UNPHASED POPSET TESTING"
	# check on the unphased statistics
	up_popset = pop_classes.UnphasedPopSet(phased_popset=my_popset)
	print "UNPHASED POPSET"
	print up_popset.__str__().replace("\n", "\n\t")

	up_popset.remove_invariant_sites()
	print "\tWith invariant sites removed:"
	print up_popset.__str__().replace("\n", "\n\t")

	#instantiate attributes to be used later
	#tester("up_popset.make_iclds()")						#up_popset.icld_ancestral[this_pop] = [list of icldA per site]
	#											#up_popset.icld_derived[this_pop] = [list of icldD per site]
	tester("up_popset.make_iehhs()")
	print "\t finished calculating iehhs", (clock() - start) / 60.0, " min.***"
	tester("up_popset.make_hamming_matrix()")				#up_popset.hamming_matrix[site, pop1, pop2]
	# print up_popset.make_matrix_text(range(2, 7), range(25, 100))
	tester("up_popset.make_pis()")						#up_popset.pi_within[this_pop] = [list of pis by site]
	tester("up_popset.make_singletons()")					#up_popset.singletons[this_pop] = [list of singletons per site]
	tester("up_popset.make_derived_counts()")				#up_popset.derived_counts[this_pop] = [list of derived counts per site]
	#calculate things important only for more than 1 pop
	if len(up_popset.pop_names) > 1:
		print "\tcalculating for multiple pops", (clock() - start) / 60.0, " min.***"
		tester("my_deltadafs = up_popset.delta_dafs()")
		tester("my_global_fst = up_popset.global_fst()")
		tester("my_fsts = up_popset.fsts(this_sim.pop_pairs)")
		tester("my_xp_iehh_Ds = up_popset.xp_iehh_Ds(this_sim.pop_pairs)")
		tester("my_xp_iehh_Ss = up_popset.xp_iehh_Ss(this_sim.pop_pairs)")
		#tester("my_xp_iclds = up_popset.xp_iclds(this_sim.pop_pairs)")
	#LD based stats
	print "\tcalculating iEHH based statistics", (clock() - start) / 60.0, " min.***"
	tester("my_ihss = up_popset.iHSs()")
	tester("my_delta_iehh_Ds = up_popset.delta_iehh_Ds()")
	#tester("my_ih_clds = up_popset.ih_clds()")
	#tester("my_delta_iclds = up_popset.delta_iclds()")
	#tester("my_cld_halfs = up_popset.cld_halfs()")

	# these dictionaries will hold all of the values that are calculated for every site in the loops below
	all_dxys = {}				# window: [list of {pop: Dxy} at each site]
	all_d12s = {}				# window: [list of {pop: D12} at each site]
	all_pis = {}				# window: [list of {pop: pi} at each site]
	all_tajima_ds = {}			# window: [list of {pop: Tajima's D} at each site]
	all_fstars = {}				# window: [list of {pop: F*} at each site]
	all_dstars = {}				# window: [list of {pop: D*} at each site]
	all_hs = {}					# window: [list of {pop: (theta_H, fay_wu_H)} at each site]
	windowed_haplogroups = {}	# window: [list of Haplogroups objects for each site]
	for w in args_window_size:
		windowed_haplogroups[w] = []	#will be a list of Haplogroups objects, one for each site
		all_dxys[w] = []				#will be a list of Dxy dictionaries, one for each site, consisting of pop:siteDxy pairs
		all_d12s[w] = []				#will be a list of D12 dictionaries, one for each site, consisting of pop:siteD12 pairs

	print "\tcalculating haplogrouped stats", (clock() - start) / 60.0, " min.***"
	#classify all of the haplogroups and get all window indexes for later use
	for (posN, dpos, bpos) in zip(up_popset.common_site_indexes, up_popset.decimal_positions[up_popset.common_site_indexes], up_popset.bp_positions[up_popset.common_site_indexes]):
		print >> hgroup_log, "//\n" + "rep:%d\tpos_index:%d\tdecimal_pos:%f\tbp_pos:%d" % (repN, posN, dpos, bpos)
		# make a popset filtered for the derived allele at this site (use for haplogrouping etc.)
		tester("derived_filter = up_popset.get_people_with_allele_at_site(posN, pop_gen.DERIVED_ALLELE, strict=True)")
		print >> hgroup_log, derived_filter
		#[/later] start with largest window and then cut it down from there, only calculating once?
		# loop over the window sizes you're using
		for (scaledW, unscaledW) in zip(SCALED_WINDOW_SIZE, args_window_size):
			#get the positions to keep for the window desired
			tester("kp, ki = up_popset.get_window_sites(dpos, scaledW)")
			#print >> hgroup_log, "kp: ", kp
			#print >> hgroup_log, "ki: ", ki
			#classify the haplogroups according to this site and this window
			tester("my_haplogroups = up_popset.get_haplogroups(ki, people_indexes=derived_filter, haplogroup_log=hgroup_log)")
			print >> hgroup_log, my_haplogroups
			# add this to the list of Haplogroups objects
			tester("windowed_haplogroups[unscaledW].append(my_haplogroups)")
			# calculate Dxy in the derived popset using the haplogrouped indices and window_indexes
			tester("all_dxys[unscaledW].append(up_popset.window_dxys(windowed_haplogroups[unscaledW][-1]))")
			# calculate D12 in the whole popset using the haplogrouped indices and window keepIs (in the hapset)
			tester("all_d12s[unscaledW].append(up_popset.window_d12s(windowed_haplogroups[unscaledW][-1]))")

	print "\tcalculating windowed statistics", (clock() - start) / 60.0, " min.***"
	# haplogroup based stats
	for w in args_window_size:
		#these list comprehensions make a list of all the sites' stat values in one window (like doing a for loop over sites)
		#allele spectrum based stats
		tester("all_pis[w] = [up_popset.window_pi(hg.window_indexes) for hg in windowed_haplogroups[w]]")
		tester("all_tajima_ds[w] = [up_popset.window_tajima_d(hg.window_indexes) for hg in windowed_haplogroups[w]]")
		tester("all_dstars[w] = [up_popset.window_d_star(hg.window_indexes) for hg in windowed_haplogroups[w]]")
		tester("all_fstars[w] = [up_popset.window_f_star(hg.window_indexes) for hg in windowed_haplogroups[w]]")
		tester("all_hs[w] = [up_popset.window_fay_wu_h(hg.window_indexes) for hg in windowed_haplogroups[w]]")

	# test the printing function
	my_stat_dict = {"global.fst" : my_global_fst,
					"fst" : my_fsts,
					"delta.daf" : my_deltadafs,
					"hgroups" : windowed_haplogroups,
					"ihs" : my_ihss,
					"delta_iehh_d" : my_delta_iehh_Ds,
					"xp_iehh_d" : my_xp_iehh_Ds,
					"xp_iehh_s" : my_xp_iehh_Ss,
					"d12" : all_d12s,
					"d.star" : all_dstars,
					"dxy" : all_dxys,
					"h.output" : all_hs,
					"f.star" : all_fstars,
					"tajima.d" : all_tajima_ds,
					"window.pi" : all_pis}
	tester("up_popset.save_stats(PREFIX, this_sim.pop_pairs, 0, my_stat_dict, args_window_size)")
