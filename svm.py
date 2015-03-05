import sys
import re
from sklearn import svm

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

def main():
	samples_file, predictions_file = get_args()
	sample_ids, sample_features, sample_labels, max_length = read_samples(samples_file)
	prediction_ids, prediction_features = read_predictions(predictions_file, max_length)
	
	classifier = svm.SVC()
	classifier.fit(sample_features, sample_labels)
	predictions = classifier.predict(prediction_features)
	
	print_predictions(prediction_ids, predictions)

def get_args():
	args = sys.argv
	num_args = len(args)
	if num_args == 3:
		return args[1:num_args]
	else:
		error('Usage: python svm.py <samples file> <predictions file>')

def error(message):
	print(message)
	exit(1)

def read_samples(samples_file):
	samples = open(samples_file, 'r').readlines()
	parsed_samples = map(parse_samples_line, samples)
	filtered_samples = filter(lambda epitope: epitope is not None, parsed_samples)
	if not filtered_samples:
		error("Couldn't find any valid samples")
	samples, features, labels = zip(*filtered_samples)
	max_length = max_sequence_length(features)
	return samples, list(features), labels, max_length

def parse_samples_line(line):
	if line.startswith('Epitope'):
		return None
	cells = line.split(',')
	if len(cells) < 55:
		error('There must be at least 55 columns per sample line, except for two lines of headers')
	epitope_id, epitope_type, sequence, label = map(strip_quotes, (cells[0], cells[10], cells[11], cells[54]))
	if re.search('.* peptide', epitope_type) is None or not sequence.isalpha():
		return None
	vectorized_sequence = vectorize_sequence(sequence)
	return [epitope_id, vectorized_sequence, label]

def strip_quotes(string):
	if string[0] == '"' and string[-1] == '"':
		return string[1:-1]
	else:
		return string

def vectorize_sequence(sequence):
	return map(lambda amino_acid: AMINO_ACIDS.index(amino_acid), sequence)

def max_sequence_length(sequences):
	return reduce(lambda length, sequence: max(length, len(sequence)), sequences, 0)

def read_predictions(predictions_file, max_length):
	predictions = parse_lines(predictions_file)
	vectorized_predictions = map(vectorize_sequence, predictions)
	prediction_ids = list(range(0, len(predictions)))
	return prediction_ids, vectorized_predictions

def parse_lines(input_file):
	lines = open(input_file, 'r').readlines()
	nonempty_lines = filter(lambda line: line, lines)
	parsed_lines = map(lambda line: line.strip(), nonempty_lines)
	return parsed_lines

def print_predictions(prediction_ids, predictions):
	num_predictions = len(prediction_ids)
	if num_predictions != len(predictions):
		error('There must be the same number of prediction IDs as predictions')
	for i in range(0, num_predictions):
		print(str(prediction_ids[i]) + ': ' + predictions[i])

main()
