import argparse
import logging
import os
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
	parser = argparse.ArgumentParser(description='Identify HVRs from a coverage file')
	parser.add_argument('--coverage', '-c', dest='coverage_file', required=True,
	                    help='The coverage file generated from calculate_coverage.py')
	parser.add_argument('--output', '-o', dest='output_folder', required=True,
	                    help='Output folder where the files will be created (will be generated if it doesn\'t exist)')
	parser.add_argument('--hvr_threshold', type=float, dest='hvr_threshold', default=0.2,
	                    help='The fraction of median coverage that defines an HVR')
	parser.add_argument('--hvr_min_length', type=int, dest='hvr_min_len', default=500,
	                    help='The minimum length of an HVR')
	parser.add_argument('--sliding_window_length', type=int, dest='sliding_window_length', default=500,
	                    help='The length of the sliding window to use')
	parser.add_argument('--sliding_window_step_size', type=int, dest='sliding_window_step', default=1,
	                    help='The length of the sliding window step to use')

	parser.add_argument('--threads', '-t', dest='threads', default=16, type=int,
	                    help='Number of threads to use')
	parser.add_argument('--log', '-l', dest='logfile', default='find_viral_HVRs.log')
	parser.add_argument('--overwrite', dest='overwrite', type=bool, default=False)
	parser.add_argument('--conda_env', dest='conda_env', default='calculate.viral.abundance')
	global args
	args = parser.parse_args()
	global logger
	logger = logging.getLogger()
	logger.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fileHandler = logging.FileHandler(args.logfile)
	fileHandler.setLevel(logging.DEBUG)
	fileHandler.setFormatter(formatter)
	consoleHandler = logging.StreamHandler()
	consoleHandler.setLevel(logging.DEBUG)
	consoleHandler.setFormatter(formatter)
	logger.addHandler(fileHandler)
	logger.addHandler(consoleHandler)

	create_output_dir()
	plt.style.use('seaborn-whitegrid')

	df = pd.read_csv(args.coverage_file, sep='\t', header=None, names=['sample', 'contig', 'loc', 'depth'])
	results = applyParallel(df.groupby(['sample', 'contig']), find_HVRs)
	#results = df.groupby(['sample', 'contig']).apply(find_HVRs)
	#results.reset_index(inplace=True, drop=True)
	final_results = pd.concat(results)
	final_results.to_csv(f'{args.output_folder}/hvrs.tsv', sep='\t', index=False)


def find_HVRs(group):
	results = []
	logging.info(f"processing {group['sample'].unique()[0]} in sample {group['contig'].unique()[0]}")
	sliding_coverage = sliding_coverage_window(group['depth'],
	                                           args.sliding_window_length,
	                                           args.hvr_threshold,
	                                           args.sliding_window_step)
	coords = identify_island_coords(sliding_coverage[2,])
	for i in range(coords.shape[0]):
		if coords[i,1] - coords[i,0] +1 >= args.hvr_min_len:
			results.append((group['sample'].unique()[0],
			                group['contig'].unique()[0],
			                coords[i,0],
			                coords[i,1],
			                coords[i,1] - coords[i,0]+1))

	df = pd.DataFrame(results, columns=['sample', 'contig', 'hvr_start', 'hvr_end', 'hvr_length'])
	df = df[df['hvr_start'] >= args.sliding_window_length]
	df = df[df['hvr_end'] <= group['depth'].size - args.sliding_window_length]
	logging.info(f"processed {group['sample'].unique()[0]} in sample {group['contig'].unique()[0]} and found {df.shape[0]} HVRs")
	return df


def applyParallel(dfGrouped, func):
	with Pool(cpu_count()) as p:
		ret_list = p.map(func, [group for name, group in dfGrouped])
	return ret_list


def create_output_dir():
	try:
		os.mkdir(args.output_folder)
		logger.debug('Created output folder at {}'.format(args.output_folder))
	except OSError:
		logger.debug('Output folder at {} already exists'.format(args.output_folder))


def sliding_coverage_window(arr, window_size, min_coverage, step=0):
    """Returns an n x 2 array where the first column is the index of the start of the window and the
    second column is an array of ones or zeros indicating the presence of a putative HVR.

    :param arr : input array.
    :type arr: numpy.ndarray
    :param window_size: size of sliding window.
    :type window_size: int
    :param min_coverage: The minimum coverage
    :type min_coverage: int

    :param step: step size of sliding window. If 0, step size is set to obtain
        non-overlapping contiguous windows (that is, step=window_size).
        Defaults to 0.
    :type step: int

    :return: array
    :rtype: numpy.ndarray
    """
    n_obs = arr.shape[0]

    # validate arguments
    if window_size > n_obs:
        raise ValueError(
            "Window size must be less than or equal "
            "the size of array in first dimension."
        )
    if step < 0:
        raise ValueError("Step must be positive.")

    n_windows = 1 + int(np.floor( (n_obs - window_size) / step))

    obs_stride = arr.strides[0]
    windowed_row_stride = obs_stride * step

    new_shape = (n_windows, window_size) + arr.shape[1:]
    new_strides = (windowed_row_stride, ) + arr.strides

    strided = np.lib.stride_tricks.as_strided(
        arr,
        shape=new_shape,
        strides=new_strides,
    )
    return np.array([np.arange(strided.shape[0]),
                     np.median(strided, axis=1),
                     np.median(strided, axis=1) < min_coverage * np.median(arr)])

def identify_island_coords(a):
    a_ext = np.concatenate(( [0], a, [0] ))
    idx = np.flatnonzero(a_ext[1:] != a_ext[:-1])
    idx = idx.reshape((int(len(idx)/2), 2)).T
    idx[1,] = idx[1,] + (int(args.sliding_window_length/2))
    idx[0,] = idx[0,] + (int(args.sliding_window_length/2))


    return idx.T

if __name__ == "__main__":
	main()
