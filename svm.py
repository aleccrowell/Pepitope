import sys
import re
from random import randint, shuffle
from sklearn import svm, preprocessing
from sklearn.metrics import explained_variance_score
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import MultiLabelBinarizer
import numpy as np

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

def main():
	train_file, test_file, num_samples, num_iterations, split_pc = get_args()

	for iteration in range(0, 5):
		training_features, training_kds = read_samples(train_file+str(iteration)+'.txt', num_samples)
		test_features, test_kds = read_samples(test_file+str(iteration)+'.txt', num_samples)
		scaler = preprocessing.StandardScaler().fit(training_kds)
		training_kds = scaler.transform(training_kds)
		test_kds = scaler.transform(test_kds)
		print len(training_features)
		classifier = svm.SVR(C=1.0, degree=5,epsilon=0.0001, gamma=0.002, kernel='rbf')
		nonamer_features = map(lambda sequence: take_nine(sequence, 0), training_features)
		lb = preprocessing.MultiLabelBinarizer().fit(nonamer_features)
		nonamer_features = lb.transform(nonamer_features)
		#nonamer_features = sp.lil_matrix(nonamer_features)
		#enc = preprocessing.OneHotEncoder().fit(nonamer_features)
		#nonamer_features = enc.transform(nonamer_features).toarray()
		print 'fitting zeroeth'
		classifier.fit(nonamer_features, training_kds)
		print 'fitted zeroeth'
		for i in range(0, num_iterations):
			for j in range(0, len(training_features)):
				combinations = all_combinations(training_features[j])
				#combinations = enc.transform(combinations).toarray()
				combinations = lb.transform(combinations)
				#combinations = sp.lil_matrix(combinations)
				predictions = classifier.predict(combinations)
				zipped = zip(combinations, predictions)
				nonamer_features[j] = filter(lambda z: z[1]==min(predictions), zipped)[0][0]
			print 'fitting '+str(i+1)
			classifier.fit(nonamer_features, training_kds)
			print 'fitted '+str(i+1)

			predictions = []
			for feature in test_features:
				test_comb = all_combinations(feature)
				#test_comb = enc.transform(test_comb).toarray()
				test_comb = lb.transform(test_comb)
				#test_comb = sp.lil_matrix(test_comb)
				feature_predictions = classifier.predict(test_comb)
				predictions.append(min(feature_predictions))
			evs = explained_variance_score(test_kds, predictions)
			mae = mean_absolute_error(test_kds, predictions)
			print evs
			print mae


def get_args():
	args = sys.argv
	num_args = len(args)
	if num_args == 6:
		train_file, test_file, num_samples, num_iterations, split_pc = args[1:6]
		if num_samples == 'all':
			num_samples = sys.maxint
		else:
			try:
				num_samples = int(num_samples)
			except:
				error("The number of samples must be an integer or 'all'")
		try:
			split_pc = int(split_pc)
		except:
			error('The percentage to split must be an integer')
		try:
			num_iterations = int(num_iterations)
		except:
			error('The number of iterations must be an integer')
		return train_file, test_file, num_samples, num_iterations, split_pc
	else:
		error('Usage: python svm.py <samples file> <predictions file> <number of samples> <number of iterations> <use qualitative?>')

def error(message):
	print(message)
	exit(1)

def read_samples(samples_file, num_samples):
	sample_lines = open(samples_file, 'r').readlines()
	selected_sample_lines = sample_lines[0:min(num_samples, len(sample_lines))]
	parsed_samples = map(parse_samples_line, selected_sample_lines)
	features, kd = map(list, zip(*parsed_samples))
	return features, kd

def parse_samples_line(line):
	if not line.startswith('human'):
		return none
	cells = line.split('\t')
	sequence, kd = map(strip_quotes, (cells[4], cells[6].strip()))
	if len(sequence) < 9:
		return none
	vectorized_sequence = vectorize_sequence(sequence)
	return [vectorized_sequence, float(kd)]

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
		print(str(prediction_ids[i]) + ': ' + str(predictions[i]))

def split_data(sample_features, sample_labels, split_pc):
	zipped = zip(sample_features, sample_labels)
	shuffle(zipped)
	divider = int(split_pc*len(zipped)/100)
	training, test = zipped[:divider], zipped[divider:]
	training_features, training_labels = zip(*training)
	test_features, test_labels = zip(*test)
	return training_features, training_labels, test_features, test_labels

main()
