import sys
import re
from random import randint
from sklearn import svm

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

def main():
	samples_file, predictions_file, num_samples, num_iterations, use_qualitative = get_args()
	sample_ids, sample_features, sample_labels, sample_kds = read_samples(samples_file, num_samples)
	prediction_ids, prediction_features = read_predictions(predictions_file)
	
	nonamer_features = map(lambda sequence: take_nine(sequence, 0), sample_features)
	classifier = svm.SVC()
	classifier.fit(nonamer_features, sample_kds)
	for i in range(0, num_iterations):
		for j in range(0, len(sample_features)):
			combinations = all_combinations(sample_features[j])
			predictions = classifier.predict(combinations)
			zipped = zip(combinations, predictions)
			best = filter(lambda z: z[1]==min(predictions), zipped)[0]
			nonamer_features[j] = best
		classifier.fit(nonamer_features, sample_kds)

		predictions = []
		for feature in prediction_features:
			feature_predictions = classifier.predict(all_combinations(feature))
			if 'Positive' in feature_predictions:
				predictions.append('Positive')
			else:
				predictions.append('Negative')
		print_predictions(prediction_ids, predictions)
		print

def get_args():
	args = sys.argv
	num_args = len(args)
	if num_args == 6:
		return args[1:num_args]
	else:
		error('Usage: python svm.py <samples file> <predictions file> <number of samples> <number of iterations> <use qualitative?>')

def error(message):
	print(message)
	exit(1)

def read_samples(samples_file, num_samples):
	sample_lines = open(samples_file, 'r').readlines()
	selected_sample_lines = sample_lines[0:min(num_samples, len(sample_lines))]
	parsed_samples = map(parse_samples_line, selected_sample_lines)
	filtered_samples = filter(lambda epitope: epitope is not None, parsed_samples)
	if not filtered_samples:
		error("Couldn't find any valid samples")
	samples, features, labels, kd = map(list, zip(*filtered_samples))
	return samples, features, labels, kd

def parse_samples_line(line):
	if line.startswith('Epitope'):
		return None
	cells = line.split(',')
	if len(cells) < 55:
		return None
	epitope_id, epitope_type, sequence, label, assay_group, units, kd = map(strip_quotes, (cells[0], cells[10], cells[11], cells[54], cells[51], cells[52], cells[55]))
	if re.search('.* peptide', epitope_type) is None or not sequence.isalpha():
		return None
	vectorized_sequence = vectorize_sequence(sequence)
	sanitized_label = sanitize_label(label)
	if vectorized_sequence is None or len(vectorized_sequence) < 9 or sanitized_label is None:
		return None
	if assay_group not in ['dissociation constant KD (~IC50)', 'dissociation constant KD (~EC50)', 'half maximal inhibitory concentration (IC50)', 'dissociation constant KD', 'Competition (or equilibrium binding) approximating KD', 'half maximal effective concentration (EC50)']:
		return None
	if units not in ['nM', 'IC50 nM', 'KD nM']:
		return None
	return [epitope_id, vectorized_sequence, sanitized_label, kd]

def strip_quotes(string):
	if string[0] == '"' and string[-1] == '"':
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
		print(str(prediction_ids[i]) + ': ' + predictions[i])

main()
