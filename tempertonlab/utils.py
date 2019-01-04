import collections
import gzip
import logging
import sys
from logging.handlers import TimedRotatingFileHandler

import itertools
import numpy as np
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from kpal.klib import Profile
from pandas import DataFrame
from sklearn.preprocessing import StandardScaler


def unzip_file(filename):
	"""
	Unzips a file and returns its unzipped name
	:param filename: The file to be unzipped
	:return: The name of the unzipped file
	"""
	subprocess.call(['unpigz', filename])
	return filename.split('.gz')[0]


def zip_file(filename):
	"""
	Unzips a file and returns its unzipped name
	:param filename: The file to be unzipped
	:return: The name of the unzipped file
	"""
	subprocess.call(['pigz', filename])
	return '%s.gz' % filename


class KmerUtility:
	"""class to manage generation of kmers

	"""

	def __init__(self, kmer_size=5):
		self.kmer_size = kmer_size

		self.kmers = self.generate_kmers()

	def generate_kmers(self):
		"""
		Generates kmers
		:return: a list of kmers
		"""
		bases = ['A', 'T', 'G', 'C']
		kmer_list = []
		for i in [''.join(p) for p in itertools.product(bases, repeat=self.kmer_size)]:
			if str(Seq(i).reverse_complement()) not in kmer_list:
				kmer_list.append(i)
		return kmer_list

	def get_kmer_profile(self, record):
		"""
		Returns a kmer profile for a record
		:param record: The sequence record
		:return: a list of kmer counts
		"""
		p = Profile.from_sequences([str(record.seq)], length=self.kmer_size, name=record.id)
		p.balance()
		counts = [p.counts[p.dna_to_binary(x)] for x in self.kmers]
		return counts

	def count_kmers(self, sequence_file, filetype='fastq', rescale=True):
		"""
		Counts the kmers per sequence for a given sequence file
		:param sequence_file: The file to be counted
		:param filetype: fasta or fastq
		:param rescale: If this is true, return values that have been standardised across the contigs using the StandardScaler
		:return: a pandas dataframe with one row per sequence and the associated sequence counts
		"""
		profiles = {}
		if filetype not in ('fasta', 'fasta.gz', 'fastq', 'fastq.gz'):
			raise ValueError()

		if filetype == 'fastq' or filetype == 'fasta':
			with open(sequence_file, 'r') as handle:
				for record in SeqIO.parse(handle, filetype):
					counts = self.get_kmer_profile(record)
					profiles[record.id] = counts

		elif filetype == 'fastq.gz':
			with gzip.open(sequence_file, 'rt') as handle:
				for record in SeqIO.parse(handle, 'fastq'):
					counts = self.get_kmer_profile(record)
					profiles[record.id] = counts
		elif filetype == 'fasta.gz':
			with gzip.open(sequence_file, 'rt') as handle:
				for record in SeqIO.parse(handle, 'fasta'):
					counts = self.get_kmer_profile(record)
					profiles[record.id] = counts

		df: DataFrame = pd.DataFrame.from_dict(profiles, columns=self.kmers, orient='index')
		if rescale:
			standard_scaler = StandardScaler()
			df = df.astype('float')
			x_scaled = standard_scaler.fit_transform(df.values)

			return pd.DataFrame(x_scaled, columns=df.columns, index=df.index)
		else:
			return df


class StatsUtility:

	def __init__(self, seed=42):
		self.seed = seed

	@staticmethod
	def ecdf(data):
		"""
		Returns ecdf data for visualising distributions
		:param data: the data to be visualised
		:return: the x and y values and a seaborn scatterplot
		"""
		n = len(data)
		x = np.sort(data)
		y = np.arange(start=1, stop=n + 1) / n
		return x, y

	def permutation_sample(self, data1, data2):
		"""Generate a permutation sample from two data sets.
		"""

		# Concatenate the data sets: data
		data = np.concatenate((data1, data2))

		# Permute the concatenated array: permuted_data
		permuted_data = np.random.permutation(data)

		# Split the permuted array into two: perm_sample_1, perm_sample_2
		perm_sample_1 = permuted_data[:len(data1)]
		perm_sample_2 = permuted_data[len(data1):]

		return perm_sample_1, perm_sample_2

	def draw_perm_reps(self, data_1, data_2, func, size=1):
		"""Generate multiple permutation replicates."""

		# Initialize array of replicates: perm_replicates
		perm_replicates = np.empty(size)

		for i in range(size):
			# Generate permutation sample
			perm_sample_1, perm_sample_2 = self.permutation_sample(data_1, data_2)

			# Compute the test statistic
			perm_replicates[i] = func(perm_sample_1, perm_sample_2)

		return perm_replicates

	def calculate_bootstrapped_median_ci(self, data, n_iterations=10000, ci=95.0):
		"""
		calculate confidence intervals around a bootstrapped median
		:param ci: The confidence interval to be calculated.
		:param data: The data to be analysed
		:param n_iterations: The number of iterations for calculating the bootstrap
		:return: A namedtuple
		"""
		medians = self.generate_n_bootstrap_replicates(data, np.median, n=n_iterations)

		ci_frac = (100 - ci) / 2
		percentiles = np.percentile(medians, [ci_frac, 100 - ci_frac])
		MedianCI = collections.namedtuple('MedianCI', 'median min max CI bootstraps')
		return MedianCI(median=np.median(medians),
		                min=percentiles[0],
		                max=percentiles[1],
		                CI=ci,
		                bootstraps=medians)

	def bootstrap_replicate_1d(self, data, func):
		"""
		Helper method to generate bootstrapped hacker stats.
		:param data: The 1d data to use
		:param func: The function to calculate the statistics
		:return: The bootstrapped replicate
		"""
		bs_sample = np.random.choice(data, len(data))
		return func(bs_sample)

	def generate_n_bootstrap_replicates(self, data, func, n=10000):
		results = np.empty(n)
		for i in range(n):
			results[i] = self.bootstrap_replicate_1d(data, func)
		return results

	def draw_bs_pairs_linreg(self, x, y, size=10000):
		"""Perform pairs bootstrap for linear regression."""

		# Set up array of indices to sample from: inds
		inds = np.arange(len(x))

		# Initialize replicates: bs_slope_reps, bs_intercept_reps
		bs_slope_reps = np.empty(size)
		bs_intercept_reps = np.empty(size)

		# Generate replicates
		for i in range(size):
			bs_inds = np.random.choice(inds, size=len(inds))
			bs_x, bs_y = x[bs_inds], y[bs_inds]
			# noinspection PyTupleAssignmentBalance
			bs_slope_reps[i], bs_intercept_reps[i] = np.polyfit(bs_x, bs_y, 1)

		return bs_slope_reps, bs_intercept_reps


class PlotUtility:
	seaborn_std_style = {'axes.axisbelow': True,
	                     'axes.edgecolor': '.8',
	                     'axes.facecolor': 'white',
	                     'axes.grid': True,
	                     'axes.labelcolor': '.15',
	                     'axes.labelsize': '14',
	                     'axes.titlesize': '16',
	                     'axes.spines.bottom': True,
	                     'axes.spines.left': True,
	                     'axes.spines.right': True,
	                     'axes.spines.top': True,
	                     'figure.facecolor': 'white',
	                     'font.family': ['Times New Roman'],
	                     'font.serif': ['Times New Roman'],
	                     'grid.color': '.8',
	                     'grid.linestyle': '-',
	                     'image.cmap': 'rocket',
	                     'lines.solid_capstyle': 'round',
	                     'patch.edgecolor': 'w',
	                     'patch.force_edgecolor': True,
	                     'text.color': '.15',
	                     'xtick.bottom': False,
	                     'xtick.color': '.15',
	                     'xtick.direction': 'out',
	                     'xtick.top': False,
	                     'ytick.color': '.15',
	                     'ytick.direction': 'out',
	                     'ytick.left': False,
	                     'ytick.right': False
	                     }


class BufferedOutputHandler:
	buffer_size = 10000

	def __init__(self, output_file, log_file, output_type='fastq'):
		self.output_file = output_file
		self.output_handle = open(output_file, 'a')
		self.output_type = output_type
		self.buffer = []
		self.total_count = 0
		self.my_logger = LoggerUtility('BufferedOutputHandler_%s' % output_file, log_file).get_logger()

		self.my_logger.debug('Opened up the buffered file')

	def add_record(self, record):
		"""
		Adds a record to the buffer. If the number of records exceeds the buffer size,
		the records are written out to save memory.
		:param record: The record to be added
		:return: None
		"""
		self.buffer.append(record)
		if len(self.buffer) >= self.buffer_size:
			count = SeqIO.write(self.buffer, self.output_handle, self.output_type)
			self.output_handle.flush()
			self.total_count += self.buffer_size
			self.buffer = []
			self.my_logger.debug('Buffering %i reads (total %i)' % (count, self.total_count))

	def close_out(self, zip_output=True):
		"""
		Closes out the handle to write out records
		:param zip_output: Whether or not to zip the output
		:return: None
		"""
		SeqIO.write(self.buffer, self.output_handle, self.output_type)
		self.output_handle.flush()
		self.output_handle.close()
		self.my_logger.debug('Closing out (%i reads added, total %i)' % (len(self.buffer),
		                                                                 self.total_count + len(self.buffer)))
		if zip_output:
			zip_file(self.output_file)


class LoggerUtility:

	def __init__(self, logger_name, log_file, level=logging.DEBUG):
		self.FORMATTER = logging.Formatter("%(asctime)s — %(name)s — %(levelname)s — %(message)s")
		self.logger = logging.getLogger(logger_name)
		self.logger.setLevel(level)  # better to have too much log than not enough
		self.logger.addHandler(self.get_console_handler())
		self.logger.addHandler(self.get_file_handler(log_file))
		self.logger.propagate = False

	def get_console_handler(self):
		"""
		Creates a logger to go to the console
		:return: console handle
		"""
		console_handler = logging.StreamHandler(sys.stdout)
		console_handler.setFormatter(self.FORMATTER)
		return console_handler

	def get_file_handler(self, log_file):
		"""
		Creates a file handler to write to a file
		:param log_file: The file to write to
		:return: file handler
		"""
		file_handler = TimedRotatingFileHandler(log_file, when='midnight')
		file_handler.setFormatter(self.FORMATTER)
		return file_handler

	def get_logger(self):
		return self.logger
