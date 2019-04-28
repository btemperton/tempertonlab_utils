from unittest import TestCase

from Bio import SeqIO

from tempertonlab_utils import GenomicsUtility


class TestRenameFastas(TestCase):

	def setUp(self):
		self.genomics_utility = GenomicsUtility()
		self.reads = [x for x in SeqIO.parse('test.fa', 'fasta')]

	def test_rename_fastas_basic(self):
		renamed = self.genomics_utility.rename_fasta(self.reads, 'new')
		self.assertEqual(3, len(renamed))
		self.assertEqual('new_0000001', renamed[0].id)
		self.assertEqual('new_0000002', renamed[1].id)
		self.assertEqual('new_0000003', renamed[2].id)

	def test_rename_fastas_min_length(self):
		renamed = self.genomics_utility.rename_fasta(self.reads, 'new', minLength=10)
		self.assertEqual(1, len(renamed))
		self.assertEqual('new_0000001', renamed[0].id)

	def test_rename_fastas_max_length(self):
		renamed = self.genomics_utility.rename_fasta(self.reads, 'new', maxLength=8)
		self.assertEqual(2, len(renamed))
		self.assertEqual('new_0000001', renamed[0].id)
		self.assertEqual('new_0000002', renamed[1].id)

	def test_rename_fastas_min_max_length(self):
		renamed = self.genomics_utility.rename_fasta(self.reads, 'new', minLength=6, maxLength=8)
		self.assertEqual(2, len(renamed))
		self.assertEqual('new_0000001', renamed[0].id)
		self.assertEqual('new_0000002', renamed[1].id)

	def test_rename_fastas_len_error_length(self):
		with self.assertRaises(Exception) as context:
			self.genomics_utility.rename_fasta(self.reads, 'new', minLength=8, maxLength=6)
		self.assertTrue('MinLength is greater than Maxlength' in str(context.exception))
