from unittest import TestCase

import pandas as pd

from tempertonlab_utils import ClusterUtility


class TestClusterGenomes(TestCase):
	def setUp(self):
		self.cluster_utility = ClusterUtility('tmp.log')
		self.expected_df = pd.read_csv('baby.cover.csv').sort_values(by=['Subject', 'Query']).reset_index(drop=True)
		self.baby_coords_file = 'baby.coords'

	def test_baby_coverage(self):
		df = self.cluster_utility.create_cover_dataframe(self.baby_coords_file)
		pd.testing.assert_frame_equal(self.expected_df,
		                              df.sort_values(by=['Subject', 'Query']).reset_index(drop=True),
		                              check_like=True)
