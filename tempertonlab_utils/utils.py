import collections
import gzip
import itertools
import logging
import subprocess
import sys
import timeit
from logging.handlers import TimedRotatingFileHandler
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from kpal.klib import Profile
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

		df = pd.DataFrame.from_dict(profiles, columns=self.kmers, orient='index')
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

	def draw_bs_reps(self, data, func, size=1):
		"""Draw bootstrap replicates."""

		# Initialize array of replicates: bs_replicates
		bs_replicates = np.empty(size)

		# Generate replicates
		for i in range(size):
			bs_replicates[i] = self.bootstrap_replicate_1d(data, func)

		return bs_replicates

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

	def draw_bs_pairs(self, x, y, func, size=1):
		"""Perform pairs bootstrap for a single statistic."""

		# Set up array of indices to sample from: inds
		inds = np.arange(len(x))

		# Initialize replicates: bs_replicates
		bs_replicates = np.empty(size)

		# Generate replicates
		for i in range(size):
			bs_inds = np.random.choice(inds, size=len(inds))
			bs_x, bs_y = x[bs_inds], y[bs_inds]
			bs_replicates[i] = func(bs_x, bs_y)

		return bs_replicates

	def shift_means(self, x, y):
		"""
		This function is used to shift the means of two datasets, so that they can
		be used to test the hypothesis that the observed difference in means between two datasets
		could happen by chance. This is different to testing whether the distributions are the same
		(which can be done by permuting the combined dataset and dividing it up into new permuted
		datasets.
		:param x: dataset
		:param y: dataset
		:return:
		"""
		mean_combined = np.mean(np.concatenate((x, y)))
		x_shifted = x - np.mean(x) + mean_combined
		y_shifted = y - np.mean(y) + mean_combined
		return x_shifted, y_shifted

	def pearson_r(x, y):
		"""Compute Pearson correlation coefficient between two arrays."""
		# Compute correlation matrix: corr_mat
		corr_mat = np.corrcoef(x, y)

		# Return entry [0,1]
		return corr_mat[0, 1]


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


class GenomicsUtility:
	def get_N50(self, lengths):
		"""
		Returns the N50 value for a list of lengths
		:param lengths: a list of lengths
		:return: The N50 value
		"""
		sorted_list = np.sort(lengths)
		rev_sorted = sorted_list[::-1]
		mask = np.cumsum(rev_sorted) <= np.sum(lengths) / 2
		n50_id = np.sum(mask)
		return rev_sorted[n50_id]

	def rename_fasta(self, fastas, prefix, minLength=0, maxLength=0, removeShrubbery=True):
		"""
		changes the names of the fasta objects and numbers them.
		Basically an anvio function without needing anvio.
		:param fastas: a set of fastas
		:param prefix: the string to call the fastas
		:param minLength: minimum length for the record to be included
		:param maxLength: maximum length for the record to be included
		:param removeShrubbery: Whether or not to strip descriptions and names
		:return: a list of renamed fasta objects
		"""

		if minLength != 0 and maxLength != 0 and minLength > maxLength:
			raise Exception('MinLength is greater than Maxlength')

		count = 1
		rtnValues = []
		for f in fastas:
			f.id = f'{prefix}_{count:07d}'
			if removeShrubbery:
				f.name = ''
				f.description = ''

			if (minLength == 0 or len(f.seq) >= minLength) and (maxLength == 0 or len(f.seq) <= maxLength):
				rtnValues.append(f)
				count += 1
		return rtnValues


class ClusterUtility:
	def __init__(self, log_file):
		self.my_logger = LoggerUtility('cluster_genomes', log_file, level=logging.INFO).get_logger()

	def get_range(self, row):
		"""
		Converts the start and end coordinates to a list of loci so we can
		capture overlapping alignments between two contigs
		:param row: The row from a query/target group
		:return: a list of loci
		"""
		start = row['s2_start']
		end = row['s2_end']
		return list(range(start, end))

	def collate_multi_alignments(self, group):
		"""
		Parses a group of query/target alignments and converts them into a single
		match and coverage value
		:param group: A group of query/target alignments
		:return: the calculated match and real coverage
		"""
		# logger.debug('in collate_multi_alignments')
		# if group.shape[0] > 2:
		#	print(f'{group.s2_name.unique()[0]}, {group.s1_name.unique()[0]}')
		match = np.sum(group.s2_len_aln * group.pct_id / 100)
		# cover = np.sum(group.s2_len_aln)
		ranges = group.apply(self.get_range, axis=1)
		combined_ranges = [i for sublist in ranges for i in sublist]
		combined_ranges_remove_overlap = set(combined_ranges)
		# I'm not quite sure yet why I have to do group.shape[0] -1 to mimic Simon's script
		real_cover = len(combined_ranges_remove_overlap) + group.shape[0]
		self.my_logger.debug(f'{group.s2_name.unique()[0]} vs. {group.s1_name.unique()[0]}')

		return ((group.s2_name.unique()[0],
		         group.s1_name.unique()[0],
		         match,
		         real_cover,
		         group.s2_len.unique()[0],
		         group.s1_len.unique()[0]))

	def parse_cover_dataframe(self, coords_df):
		pass

	def applyParallel(self, dfGrouped, func):
		with Pool(cpu_count()) as p:
			ret_list = p.map(func, [group for name, group in dfGrouped])
		return ret_list

	def parse_coords_file(self, coords_file, multithreaded=True):
		self.my_logger.debug(f'Parsing file: {coords_file}')
		tic = timeit.default_timer()
		coords_df = pd.read_csv(coords_file, sep='\t',
		                        header=None,
		                        names=['s1_start', 's1_end',
		                               's2_start', 's2_end',
		                               's1_len_aln', 's2_len_aln',
		                               'pct_id', 's1_len', 's2_len',
		                               'cov_1', 'cov_2',
		                               's1_name', 's2_name'])

		coords_df = coords_df[coords_df.s1_name != coords_df.s2_name]
		toc = timeit.default_timer()
		self.my_logger.debug(f'Loading the coords file took {(toc - tic):.2f} seconds')
		if multithreaded:
			results = self.applyParallel(coords_df.groupby(['s1_name', 's2_name']),
			                             self.collate_multi_alignments)
		else:
			results = coords_df.groupby(['s1_name', 's2_name']).apply(self.collate_multi_alignments).tolist()

		df = pd.DataFrame(results, columns=['Subject', 'Query', 'matches', 'hit_len', 'Subject_length', 'Query_length'])
		return df

	def create_cover_dataframe(self, coords_file, multithreaded=True):
		df = self.parse_coords_file(coords_file, multithreaded=multithreaded)
		return df
