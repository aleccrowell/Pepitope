import sys
import re
import matplotlib.pyplot as plt
import numpy as np

COLORS = ['r', 'g', 'b', 'c', 'm']

def main():
	data_file, num_folds, image_file = get_args()
	cells = None
	try:
		cells = map(lambda line: filter(lambda cell: cell, re.split(' |\n', line)), open(data_file, 'r').readlines())
	except:
		error("Couldn't read '" + data_file + "'")
	num_rows = len(cells)
	try:
		for i in range(num_rows):
			for j in range(0, 4):
				cells[i][j] = int(cells[i][j])
			cells[i][4] = float(cells[i][4])
	except:
		error("Couldn't parse '" + data_file + "'")
	ax = plt.subplot(111)
	num_iterations = int(num_rows / num_folds)
	average_aucs = []
	for i in np.arange(0, num_folds*num_iterations, num_iterations):
		fold = cells[i:i+num_iterations]
		aucs = map(lambda cells: cells[4], fold)
		ax.plot(range(num_iterations), aucs, COLORS[int(i/num_iterations)]+'o-')
	ax.set_title('AUC per Iteration')
	ax.set_xlabel('Iteration')
	ax.set_ylabel('AUC score')
	plt.savefig(image_file)

def get_args():
	args = sys.argv
	num_args = len(args)
	if num_args == 4:
		data_file, num_folds, image_file = args[1], args[2], args[3]
		try:
			num_folds = int(num_folds)
		except:
			error("Couldn't parse the given number of folds, '" + num_folds + "', as an integer")
		return data_file, num_folds, image_file
	else:
		error('Usage: python iterations.py <data file to read from> <number of cross validation num_folds used> <image file to write to>')

def error(message):
	print(message)
	exit(1)

main()
