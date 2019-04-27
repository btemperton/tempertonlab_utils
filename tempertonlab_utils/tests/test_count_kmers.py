from unittest import TestCase

from scipy import stats

from tempertonlab_utils import KmerUtility


class TestCountKmers(TestCase):
	def setUp(self):
		self.kmer_utility = KmerUtility()
		self.test_counts = stats.zscore([2, 2, 4])


class TestGenerateKmers(TestCountKmers):
	def test_generate_kmers(self):
		self.assertEqual(512, len(self.kmer_utility.kmers))

	def test_count_kmers_fasta(self):
		kmers_df = self.kmer_utility.count_kmers('test.fa', filetype='fasta')

		self.assertEqual(3, kmers_df.shape[0])
		self.assertEqual(self.test_counts[0], kmers_df.loc['test_1']['AAAAA'])
		self.assertEqual(self.test_counts[2], kmers_df.loc['test_3']['AAAAA'])

	def test_count_kmers_fasta_gz(self):
		kmers_df = self.kmer_utility.count_kmers('test.fa.gz', filetype='fasta.gz')
		self.assertEqual(3, kmers_df.shape[0])
		self.assertEqual(self.test_counts[0], kmers_df.loc['test_1']['AAAAA'])
		self.assertEqual(self.test_counts[2], kmers_df.loc['test_3']['AAAAA'])

	def test_count_kmers_fastq(self):
		kmers_df = self.kmer_utility.count_kmers('test.fq', filetype='fastq')
		self.assertEqual(3, kmers_df.shape[0])
		self.assertEqual(self.test_counts[0], kmers_df.loc['test_1']['AAAAA'])
		self.assertEqual(self.test_counts[2], kmers_df.loc['test_3']['AAAAA'])

	def test_count_kmers_fastq_gz(self):
		kmers_df = self.kmer_utility.count_kmers('test.fq.gz', filetype='fastq.gz')
		self.assertEqual(3, kmers_df.shape[0])
		self.assertEqual(self.test_counts[0], kmers_df.loc['test_1']['AAAAA'])
		self.assertEqual(self.test_counts[2], kmers_df.loc['test_3']['AAAAA'])
