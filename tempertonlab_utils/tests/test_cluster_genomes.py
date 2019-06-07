import timeit
from unittest import TestCase

import pandas as pd

from tempertonlab_utils import ClusterUtility


class TestClusterGenomes(TestCase):
	def setUp(self):
		self.cluster_utility = ClusterUtility('tmp.log')
		self.expected_df = pd.read_csv('baby.cover.csv').sort_values(by=['Subject', 'Query']).reset_index(drop=True)
		self.baby_coords_file = 'baby.coords'
		self.large_coords = 'ava.nucmer.coords'

	def test_baby_coverage_single_thread(self):
		df = self.cluster_utility.create_cover_dataframe(self.baby_coords_file, multithreaded=False)
		pd.testing.assert_frame_equal(self.expected_df,
		                              df.sort_values(by=['Subject', 'Query']).reset_index(drop=True),
		                              check_like=True)

	def test_baby_coverage_multi_thread(self):
		df = self.cluster_utility.create_cover_dataframe(self.baby_coords_file, multithreaded=True)
		pd.testing.assert_frame_equal(self.expected_df,
		                              df.sort_values(by=['Subject', 'Query']).reset_index(drop=True),
		                              check_like=True)

	def test_large_coverage_multithreaded(self):
		tic = timeit.default_timer()
		df = self.cluster_utility.create_cover_dataframe(self.large_coords, multithreaded=True)
		toc = timeit.default_timer()
		print(f'That took {toc - tic} seconds for a multithreaded approach')
