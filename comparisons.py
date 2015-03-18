import sys
import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

COLORS = ['r', 'g', 'b', 'c', 'm', 'y']

def main():
	our_data_files, their_data_file, individual_plot_file, combined_plot_file = get_args()
	
	cells = None
	aucs = {}
	for data_file in our_data_files:
		match = re.search('output-(.+)', data_file)
		if match is None:
			error("Couldn't parse allele from filename")
		else:
			allele = match.group(1)
			aucs[allele] = []
		try:
			file_aucs = map(lambda line: filter(lambda cell: cell, re.split(' ', line))[-1], open(data_file, 'r').readlines())
		except:
			error("Couldn't read '" + data_file + "'")
		try:
			auc_floats = map(lambda auc: float(auc), file_aucs)
			average_auc = sum(auc_floats) / len(file_aucs)
			aucs[allele].append(average_auc)
		except:
			error("Couldn't parse '" + data_file + "'")
	
	headers = ['Pepitope SVM']
	try:
		cells = map(lambda line: filter(lambda cell: cell, re.split('\t|\n', line)), open(their_data_file, 'r').readlines())
		headers += cells[0][1:]
		for row in cells[1:]:
			aucs[row[0]] += map(float, row[1:])
	except:
		error("Couldn't parse '" + their_data_file + "'")
	
	ax = plt.subplot(111)
	patches = []
	bar_width = 0.12
	i = 0
	for allele in aucs:
		allele_aucs = aucs[allele]
		lefts = np.arange(len(allele_aucs)) + bar_width * i
		ax.bar(lefts, allele_aucs, bar_width, color=COLORS[i], label=headers[i])
		patches.append(mpatches.Patch(color=COLORS[i], label=headers[i]))
		i += 1
	ax.set_title('AUC Comparisons')
	ax.set_xlabel('Allele')
	ax.set_ylabel('AUC score')
	ax.set_yticks(np.arange(0.0, 1.1, 0.1))
	plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1))
	fig = plt.gcf()
	fig.set_size_inches(12, 10)
	plt.savefig(individual_plot_file)
	plt.clf()
	
	aucs_by_method = []
	for allele in aucs:
		allele_aucs = aucs[allele]
		for i in range(len(allele_aucs)):
			if len(aucs_by_method) <= i:
				aucs_by_method.append([])
			aucs_by_method[i].append(allele_aucs[i])
		
	ax = plt.subplot(111)
	ax.boxplot(aucs_by_method)
	ax.set_title('AUC Comparisons')
	ax.set_xlabel('Method')
	ax.set_ylabel('AUC score')
	ax.set_xticklabels(headers, rotation=45)
	plt.savefig(combined_plot_file)

def get_args():
	args = sys.argv
	num_args = len(args)
	if num_args >= 5:
		our_data_files, their_data_file, individual_plot_file, combined_plot_file = args[1:-3], args[-3], args[-2], args[-1]
		return our_data_files, their_data_file, individual_plot_file, combined_plot_file
	else:
		error('Usage: python iterations.py <data files from an SVM run> <data file summarizing other results> <image file for individual comparisons> <image file for combined comparisons>')

def error(message):
	print(message)
	exit(1)

main()
