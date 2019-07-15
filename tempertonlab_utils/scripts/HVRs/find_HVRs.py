import argparse
import logging
import os
import re
import subprocess
import warnings
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
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
	parser.add_argument('--fasta', required=True, dest='fasta_file')

	parser.add_argument('--threads', '-t', dest='threads', default=16, type=int,
	                    help='Number of threads to use')
	parser.add_argument('--log', '-l', dest='logfile', default='find_viral_HVRs.log')
	parser.add_argument('--overwrite', dest='overwrite', type=bool, default=False)
	parser.add_argument('--conda_env', dest='conda_env', default='calculate.viral.abundance')
	parser.add_argument('--pvog', dest='pvog_db', default='../../dbs/pvog.hmm', help='The location of the pVOG HMM database for prokka to use')
	global args
	args = parser.parse_args()
	global logger
	logger = logging.getLogger()

	global traces
	traces = dict()
	logger.setLevel(logging.INFO)
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fileHandler = logging.FileHandler(args.logfile)
	fileHandler.setLevel(logging.INFO)
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

	final_results = pd.concat(results)
	final_results.to_csv(f'{args.output_folder}/hvrs.tsv', sep='\t', index=False)

	anotate_genomes(args.fasta_file)

def anotate_genomes(filename):
	logger.info('Annotating genomes')

	genes_df = get_genes_with_prodigal(args.fasta_file, f'{args.output_folder}/prots.fa')

	ribosomal_df = get_ribosomal(args.fasta_file, f'{args.output_folder}/ribosomal.gff')

	hmmscan_results_df = runHMMScan(f'{args.output_folder}/prots.fa', f'{args.output_folder}/pvog.tbl')


def get_ribosomal(infile, outfile):
	if os.path.isfile(outfile):
		logger.info(f'ribosomal file {outfile} already exists')
	else:
		logger.info('Calling barrnap')
		cmd = f'barrnap --threads {args.threads} {infile}'
		stdout, stderr = execute(cmd)
		logger.info(stderr.decode('utf-8'))
		logger.debug(stdout.decode('utf-8'))
		with open(outfile, 'w') as handle:
			handle.write(stdout.decode('utf-8'))



def get_genes_with_prodigal(infile, outfile):
	if os.path.isfile(outfile):
		logger.info(f'Protein file {outfile} already exists')
	else:
		logger.info('Calling genes with prodigal')
		cmd = f"""prodigal -a {outfile} -f gff -i {infile} -m -p meta -q"""
		stdout, stderr = execute(cmd)
		logger.info(stderr.decode('utf-8'))
		logger.debug(stdout.decode('utf-8'))

	if os.path.isfile(f'{args.output_folder}/genes_df.tsv'):
		gene_df = pd.read_csv(f'{args.output_folder}/genes_df.tsv', sep='\t')
	else:
		rgx = re.compile('((\S+)_(\S+))\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+([0-9\-]+).*')
		protein_results = []
		cleaned_proteins = []
		with open(outfile, 'r') as handle:
			for record in SeqIO.parse(handle, 'fasta'):
				m = rgx.search(record.description)
				if m:
					protein_results.append((m.group(2), m.group(1), m.group(4), m.group(5), m.group(6), 'prodigal'))
				record.name = ''
				record.description = ''
				record.seq = record.seq[:-1]
				cleaned_proteins.append(record)
		SeqIO.write(cleaned_proteins, outfile, 'fasta')
		gene_df = pd.DataFrame(protein_results, columns=['contig', 'protein_id', 'start', 'end', 'frame', 'method'])
		gene_df.to_csv(f'{args.output_folder}/genes_df.tsv', sep='\t', index=False)
	return gene_df





def runHMMScan(infile, outfile):
	if os.path.isfile(outfile):
		logger.info(f'PVOG output file {args.output_folder}pvog.tbl already exists')
	else:
		logger.info('Running HMM search against pVOG')
		cmd = f'hmmsearch --tblout {outfile} --notextw -E 1e-5 --cpu {args.threads} {args.pvog_db} {infile}'
		stdout, stderr = execute(cmd)
		logger.info(stderr.decode('utf-8'))
		logger.debug(stdout.decode('utf-8'))

	hmm_results = pd.read_csv(outfile, sep='\s+', skiprows=3, comment='#',usecols=[0,2])
	return hmm_results






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

	store_trace(contig_name, sample_name, sliding_coverage)



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

def store_trace(contig_name, sample_name, sliding_coverage):
	try:
		traces[contig_name][sample_name] = sliding_coverage[1,]
	except:
		traces[contig_name] = dict()
		traces[contig_name]['coords'] = sliding_coverage[0,]
		traces[contig_name][sample_name] = sliding_coverage[1,]


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
    return idx.T

def execute(command):
	logger.info('\nExecuting {}'.format(command))
	process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
	                           shell=True)
	(stdout, stderr) = process.communicate()

	return stdout, stderr


if __name__ == "__main__":
	main()
