import sys
import re
import matplotlib.pyplot as plt

def main():
	data_files, image_file = get_args()
	data = []
	for data_file in data_files:
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
		skewness = cells[0][0] / float(cells[0][0] + cells[0][1])
		average_auc = sum(map(lambda c: c[4], cells)) / num_rows
		data.append((data_file, skewness, average_auc))
	sorted_data = sorted(data, key=lambda d: d[1])
	filenames, skewnesses, average_aucs = zip(*sorted_data)
	ax = plt.subplot(111)
	ax.plot(skewnesses, average_aucs, 'bo-')
	ax.set_title('AUC vs Skewness')
	ax.set_xlabel('Skewness (positive / total)')
	ax.set_ylabel('AUC score')
	plt.savefig(image_file)

def get_args():
	args = sys.argv
	num_args = len(args)
	if num_args >= 3:
		return args[1:-1], args[-1]
	else:
		error('Usage: python skew_v_auc.py <data files to read from> <image file to write to>')

def error(message):
	print(message)
	exit(1)

main()
