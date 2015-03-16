import sys
import re
from random import randint, shuffle
from sklearn import svm
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import MultiLabelBinarizer
import numpy as np

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

def main():
	samples_file, predictions_file, num_samples, num_iterations, split_pc = get_args()
	sample_ids, sample_features, sample_labels = read_samples(samples_file, num_samples)
	
	mlb = MultiLabelBinarizer()
	for fold in range(0, int(100/(100-split_pc))):
		training_features, training_labels, test_features, test_labels = split_data(sample_features, sample_labels, split_pc)
		if has_both(training_labels) and has_both(test_labels):
			training_trues, training_falses = len(filter(lambda l: l, training_labels)), len(filter(lambda l: not l, training_labels))
			test_trues, test_falses = len(filter(lambda l: l, test_labels)), len(filter(lambda l: not l, test_labels))
			nonamer_features = map(lambda sequence: take_nine(sequence, 0), training_features)
			classifier = svm.SVC(probability=True)
			classifier.fit(nonamer_features, training_labels)
			for i in range(0, num_iterations):
				for j in range(0, len(training_features)):
					combinations = all_combinations(training_features[j])
					predictions = classifier.predict(combinations)
					zipped = zip(combinations, predictions)
					best = filter(lambda z: z[1]==min(predictions), zipped)[0][0]
					nonamer_features[j] = best
				classifier.fit(nonamer_features, training_labels)
	
				predictions = []
				for feature in test_features:
					feature_combinations = all_combinations(feature)
					feature_predictions = classifier.predict(feature_combinations)
					if any(feature_predictions):
						prediction = feature_combinations[feature_predictions.tolist().index(True)]
					else:
						prediction = feature_combinations[0]
					predictions.append(prediction)
	
				probabilities = classifier.predict_proba(predictions)
				test_labels_mlb = mlb.fit_transform(map(lambda l: [1] if l else [0], test_labels))
				print str(training_trues)+' '+str(training_falses)+' '+str(test_trues)+' '+str(test_falses)+' '+str(roc_auc_score(test_labels_mlb, probabilities, average='micro'))

def get_args():
	args = sys.argv
	num_args = len(args)
	if num_args == 6:
		samples_file, predictions_file, num_samples, num_iterations, split_pc = args[1:6]
		if num_samples == 'all':
			num_samples = sys.maxint
		else:
			try:
				num_samples = int(num_samples)
			except:
				error("The number of samples must be an integer or 'all'")
		try:
			split_pc = int(split_pc)
			if split_pc <= 0 or split_pc >= 100:
				error('The percentage to split must be an integer greater than 0 and less than 100')
		except:
			error('The percentage to split must be an integer')
		try:
			num_iterations = int(num_iterations)
		except:
			error('The number of iterations must be an integer')
		return samples_file, predictions_file, num_samples, num_iterations, split_pc
	else:
		error('Usage: python svm.py <samples file> <predictions file> <number of samples> <number of iterations> <percent of data to train>')

def error(message):
	print(message)
	exit(1)

def read_samples(samples_file, num_samples):
	sample_lines = open(samples_file, 'r').readlines()
	selected_sample_lines = sample_lines[0:min(num_samples, len(sample_lines))]
	if samples_file.endswith(".csv"):
		parsed_samples = map(parse_iedb_line, selected_sample_lines)
	elif samples_file.endswith(".tab"):
		parsed_samples = map(parse_rta_line, selected_sample_lines)
	filtered_samples = filter(lambda epitope: epitope is not None, parsed_samples)
	if not filtered_samples:
		error("Couldn't find any valid samples")
	samples, features, kd = map(list, zip(*filtered_samples))
	return samples, features, kd

def parse_iedb_line(line):
	if line.startswith('Epitope') or line.startswith('MHC ligand ID') or line.startswith('Reference ID'):
		return None
	cells = line.split(',')
	if len(cells) < 55 or not is_floatable(cells[54]):
		return None
	epitope_id, epitope_type, sequence, label, assay_group, units, kd = map(strip_quotes, (cells[0], cells[9], cells[10], cells[53], cells[50], cells[51], cells[54]))
	if re.search('.* peptide', epitope_type) is None or not sequence.isalpha():
		return None
	vectorized_sequence = vectorize_sequence(sequence)
	sanitized_label = sanitize_label(label) # not currently used
	if vectorized_sequence is None or len(vectorized_sequence) < 9 or sanitized_label is None:
		return None
	if assay_group not in ['dissociation constant KD (~IC50)', 'dissociation constant KD (~EC50)', 'half maximal inhibitory concentration (IC50)', 'dissociation constant KD', 'Competition (or equilibrium binding) approximating KD', 'half maximal effective concentration (EC50)']:
		return None
	if units not in ['nM', 'IC50 nM', 'KD nM']:
		return None
	return [epitope_id, vectorized_sequence, float(kd)]

def is_floatable(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

def strip_quotes(string):
	if type(string) == 'str' and string[0] == '"' and string[-1] == '"':
		return string[1:-1]
	else:
		return string

def vectorize_sequence(sequence):
	if len(filter(lambda amino_acid: amino_acid in AMINO_ACIDS, sequence)) == len(sequence):
		return map(lambda amino_acid: AMINO_ACIDS.index(amino_acid), sequence)
	else:
		return None

def sanitize_label(label):
	if label in ['Positive', 'Positive-Intermediate', 'Positive-High', 'Positive-Low']:
		return 'Positive'
	elif label in ['Negative']:
		return 'Negative'
	else:
		return None

def parse_rta_line(line):
	cells = re.split('\t|,', line)
	if len(cells) >= 5:
		epitope_id = cells[0]
		vectorized_sequence = vectorize_sequence(cells[3])
		if vectorized_sequence is None or len(vectorized_sequence) < 9:
			return None
		if cells[4] == 'Negative':
			kd = False
		else:
			kd = True
		return [epitope_id, vectorized_sequence, kd]
	else:
		return None

def read_predictions(predictions_file):
	predictions = parse_lines(predictions_file)
	vectorized_predictions = map(vectorize_sequence, predictions)
	prediction_ids = list(range(0, len(predictions)))
	return prediction_ids, vectorized_predictions

def parse_lines(input_file):
	lines = open(input_file, 'r').readlines()
	nonempty_lines = filter(lambda line: line, lines)
	parsed_lines = map(lambda line: line.strip(), nonempty_lines)
	return parsed_lines

def has_both(labels):
	return any(labels) and any(map(lambda l: not l, labels))

def take_nine(sequence, start):
	return sequence[start:start+9]

def all_combinations(sequence):
	return map(lambda i: sequence[i:i+9], range(0, len(sequence)-8))

def random_prediction(pairs):
	return pairs[randint(0, len(pairs)-1)][0]

def print_predictions(prediction_ids, predictions):
	num_predictions = len(prediction_ids)
	if num_predictions != len(predictions):
		error('There must be the same number of prediction IDs as predictions')
	for i in range(0, num_predictions):
		print(str(prediction_ids[i]) + ': ' + str(predictions[i]))

def split_data(sample_features, sample_labels, split_pc):
	zipped = zip(sample_features, sample_labels)
	shuffle(zipped)
	divider = int(split_pc*len(zipped)/100)
	training, test = zipped[:divider], zipped[divider:]
	training_features, training_labels = map(list, zip(*training))
	test_features, test_labels = map(list, zip(*test))
	return training_features, training_labels, test_features, test_labels

main()
