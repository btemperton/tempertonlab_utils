import argparse
import logging
import os
import warnings
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection

warnings.simplefilter(action='ignore', category=FutureWarning)


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

	create_output_dir(args.output_folder)
	plt.style.use('seaborn-ticks')

	df = pd.read_csv(args.coverage_file, sep='\t', header=None, names=['sample', 'contig', 'loc', 'depth'])
	results = applyParallel(df.groupby(['sample', 'contig']), find_HVRs)
	#results = df.groupby(['sample', 'contig']).apply(find_HVRs)
	#results.reset_index(inplace=True, drop=True)

	final_results = pd.concat(results)
	final_results.to_csv(f'{args.output_folder}/hvrs.tsv', sep='\t', index=False)


def plot_coverage(sliding_coverage, name, sample, hvr_df):
	fig, ax = plt.subplots(figsize=(10, 4))
	ax.fill_between(sliding_coverage[0,], np.sqrt(sliding_coverage[1,]), color='green')

	boxes = [get_hvr_rect(start, width,
	                       np.max(sliding_coverage[1,]))
	          for start, width in zip(hvr_df['hvr_start'], hvr_df['hvr_length'])]
	pc = PatchCollection(boxes,facecolor='lightgrey', alpha=0.6)
	ax.add_collection(pc)
	ax.set(title=f"Coverage of {name} in {sample}",
	       ylabel=r'$\sqrt{Coverage}$',
	       xlabel='Locus (bp)')



	fig.tight_layout()
	create_output_dir(f'{args.output_folder}/{name}')
	plt.savefig(f'{args.output_folder}/{name}/{sample}.png', dpi=300)

def get_hvr_rect(start, width, height):
	rectangle = plt.Rectangle((start, 0), width, height)
	return rectangle


def find_HVRs(group):
	results = []
	sliding_coverage = sliding_coverage_window(group['depth'],
	                                           args.sliding_window_length,
	                                           args.hvr_threshold,
	                                           args.sliding_window_step)

	contig_name = group['contig'].unique()[0]
	sample_name = group['sample'].unique()[0]




	coords = identify_island_coords(sliding_coverage[2,])
	for i in range(coords.shape[0]):
		if coords[i,1] - coords[i,0] +1 >= args.hvr_min_len:
			results.append((sample_name,
			                contig_name,
			                coords[i,0],
			                coords[i,1],
			                coords[i,1] - coords[i,0]+1))

	df = pd.DataFrame(results, columns=['sample', 'contig', 'hvr_start', 'hvr_end', 'hvr_length'])
	df = df[df['hvr_start'] >= args.sliding_window_length]
	df = df[df['hvr_end'] <= group['depth'].size - args.sliding_window_length]

	plot_coverage(sliding_coverage, contig_name, sample_name, df)

	logging.info(f"processed {contig_name} in sample {sample_name} and found {df.shape[0]} HVRs")
	return df


def applyParallel(dfGrouped, func):
	with Pool(cpu_count()) as p:
		ret_list = p.map(func, [group for name, group in dfGrouped])
	return ret_list


def create_output_dir(directory_name):
	try:
		os.mkdir(directory_name)
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
    #We need to do this bit because we're taking medians.
    #idx[1,] = idx[1,] + (int(args.sliding_window_length/2))
    #idx[0,] = idx[0,] - (int(args.sliding_window_length/2))


    return idx.T



if __name__ == "__main__":
	main()
