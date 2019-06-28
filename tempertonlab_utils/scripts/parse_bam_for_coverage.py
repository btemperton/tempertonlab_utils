import argparse
import glob
import logging
import os
import subprocess
from io import StringIO
import pandas as pd
import numpy as np
from multiprocessing import Pool


def main():
	parser = argparse.ArgumentParser(description='Create mapping files for viral genomes against a set of metagenomes')
	parser.add_argument('--dir', '-d', dest='bamm_dir', required=True, help='Directory containing the bamm filtered bam files')
	parser.add_argument('--pct_cutoff', '-c', dest='pct_cutoff', default =75,
						help='The cutoff used to determine if a virus is in a sample (e.g. 75% of the genome has to be covered')
	parser.add_argument('--threads', '-t', dest='thread_count', default=16, type=int,
						help='The number of threads to use')
	parser.add_argument('--log', '-l', dest='logfile', default='parse_bam_for_coverage.log')
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

	bam_files = glob.glob(f'{os.path.abspath(args.bamm_dir)}/*.bam')

	with Pool(processes=args.thread_count) as pool:
		dfs = pool.map(parse_bams, bam_files)

	final_df = pd.concat(dfs)
	final_df.to_include = final_df.to_include.astype(int)


	final_df.to_csv(f'{os.path.abspath(args.bamm_dir)}/contig.{args.pct_cutoff}pc.details.tsv', sep='\t', index=False)



def calculate_tpmean(s):
	min, max = np.quantile(s, [0.1, 0.9])
	masked_arr = np.where((s <= min) | (s >= max), 0, s)
	logger.debug(f'Before tpmean, sum is {np.sum(s)}; After it is {np.sum(masked_arr)}')
	return np.mean(masked_arr)


def summarise_group(x):
	names = {
		'num_bases_covered': x['is_covered'].sum(),
		'total_bases': x['is_covered'].count(),
		'tpmean': calculate_tpmean(x['coverage'].to_numpy())
	}
	return pd.Series(names, index=['num_bases_covered', 'total_bases', 'tpmean'])


def parse_bams(bam_file):
	logger.debug(f'Parsing file {bam_file}')

	cmd = f'samtools depth -aa {bam_file}'
	stdout, stderr = execute(cmd)

	df = pd.read_csv(StringIO(stdout.decode('utf-8')), sep='\t', header=None, names=['contig', 'pos', 'coverage'])

	df['is_covered'] = df.coverage.astype(bool).astype(int)

	summary_df = df.groupby(['contig']).apply(summarise_group)
	summary_df.reset_index(level=0, inplace=True)
	summary_df['sample'] = os.path.basename(bam_file).strip('.bam')
	summary_df['frac_covered'] = summary_df['num_bases_covered'] / summary_df['total_bases'] *100
	summary_df['to_include'] = summary_df['frac_covered'] >= args.pct_cutoff

	logger.debug(f'''In total, for sample {os.path.basename(bam_file).strip(".bam")},
	{summary_df[summary_df.to_include].shape[0]} contigs will be included''')
	return summary_df[['sample', 'contig', 'frac_covered', 'to_include', 'tpmean']]


def execute(command):
	logger.info('\nExecuting {}'.format(command))
	process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(stdout, stderr) = process.communicate()

	return stdout, stderr


if __name__ == "__main__":
	main()
