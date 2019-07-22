import argparse
import logging
import os
import re
import subprocess
import warnings
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.SeqRecord import SeqRecord

warnings.simplefilter(action='ignore', category=FutureWarning)


def main():
	parser = argparse.ArgumentParser(description='Annotate genomes')
	parser.add_argument('--output', '-o', dest='output_folder', required=True,
	                    help='Output folder where the files will be created (will be generated if it doesn\'t exist)')
	parser.add_argument('--fasta', required=True, dest='fasta_file')

	parser.add_argument('--threads', '-t', dest='threads', default=16, type=int,
	                    help='Number of threads to use')
	parser.add_argument('--log', '-l', dest='logfile', default='find_viral_HVRs.log')
	parser.add_argument('--overwrite', dest='overwrite', type=bool, default=False)
	parser.add_argument('--conda_env', dest='conda_env', default='calculate.viral.abundance')
	parser.add_argument('--pvog', dest='pvog_db', default='/Users/bt273/Google Drive/TempertonLab/tempertonlab_utils/tempertonlab_utils/dbs/pvog.hmm',
	                    help='The location of the pVOG HMM database for prokka to use')
	global args
	args = parser.parse_args()
	global logger
	logger = logging.getLogger()

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
	annotations_df = annotate_genomes(args.fasta_file)


def annotate_genomes(filename):
	logger.info('Annotating genomes')

	prodigal_genes_df = get_genes_with_prodigal(filename, f'{args.output_folder}/prodigal.tsv')

	mga_genes_df = get_genes_with_mga(filename, f'{args.output_folder}/mga.tsv')

	ribosomal_df = get_ribosomal(filename, f'{args.output_folder}/ribosomal.gff')

	#final_df = prodigal_genes_df.merge(hmmscan_results_df, how='left', on='protein_id')
	#return final_df


def create_output_dir(directory_name):
	try:
		os.mkdir(directory_name)
		logger.debug('Created output folder at {}'.format(args.output_folder))
	except OSError:
		logger.debug('Output folder at {} already exists'.format(args.output_folder))


def execute(command):
	logger.info('\nExecuting {}'.format(command))
	process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
	                           shell=True)
	(stdout, stderr) = process.communicate()

	return stdout, stderr


def get_genes_with_mga(input_file, output_file):
	if os.path.isfile(output_file):
		logger.info(f'MGA output file {output_file} already exists')
		merged_df = pd.read_csv(output_file, sep='\t')
	else:
		cmd = f'mga {input_file} -m'
		stdout, stderr = execute(cmd)
		df = parse_mga_output(stdout.decode('utf-8'))
		df['method'] = 'mga'

		contigs = SeqIO.index(input_file, 'fasta')
		predicted_proteins = [get_predicted_protein(contigs[contig], name, start, end, frame)
		                      for contig, name, start, end, frame in zip(df['contig'],
		                                                                 df['gene_id'],
		                                                                 df['start'],
		                                                                 df['end'],
		                                                                 df['frame'])]
		predicted_proteins = [x for x in predicted_proteins if x]
		SeqIO.write(predicted_proteins, f'{args.output_folder}/mga.fa', 'fasta')

		hmm_results = runHMMScan(f'{args.output_folder}/mga.fa')
		hmm_results['contig'], hmm_results['gene_id'] = hmm_results['protein_id'].str.split('_', 1).str
		hmm_results = hmm_results[['contig', 'gene_id', 'pVOG']]
		merged_df = df.merge(hmm_results, how='left', on=['contig', 'gene_id'])
		logging.info(f'Writing file {os.path.abspath(output_file)}')
		merged_df.to_csv(os.path.abspath(output_file), sep='\t', index=False)

	return merged_df


def get_predicted_protein(contig, name, start, finish, frame):
	logger.debug(f'Getting protein from {contig} at {start}-{finish}')
	subseq = contig.seq[int(start) -1: int(finish)]
	if frame is -1:
		subseq = subseq.reverse_complement()
	try:
		aa_seq = subseq.translate(table="Bacterial", cds=True)
		return SeqRecord(aa_seq, id=f'{contig.id}_{name}', description = '')
	except TranslationError as e:
		logger.warning(f'{contig.id} {start}-{finish}: {str(e)}')


def parse_mga_output(output):
	results = []
	contig = ''
	contig_rgx = re.compile('# (\S+)')
	frame_map = {'+': 1, '-': -1}
	lines = iter(output.split('\n'))
	for line in lines:
		if line.startswith('#'):
			m = contig_rgx.search(line)
			if m:
				contig = m.group(1)
				next(lines, None)
				next(lines, None)
				continue
		bits = line.split('\t')
		try:
			results.append((contig, bits[0], bits[1], bits[2], frame_map[bits[3]]))
		except IndexError:
			pass

	df = pd.DataFrame(results, columns=['contig', 'gene_id', 'start', 'end', 'frame'])

	return df


def applyParallel(dfGrouped, func):
	with Pool(cpu_count()) as p:
		ret_list = p.map(func, [group for name, group in dfGrouped])
	return ret_list


def runHMMScan(infile):

	logger.info('Running HMM search against pVOG')
	cmd = f'''hmmsearch --tblout tmp.tbl \
	--notextw \
	-E 1e-5 \
	--cpu {args.threads} \
	"{os.path.abspath(args.pvog_db)}" \
	"{os.path.abspath(infile)}"'''
	stdout, stderr = execute(cmd)
	if stderr:
		logger.info(stderr.decode('utf-8'))
	if stdout:
		logger.debug(stdout.decode('utf-8'))

	hmm_results = pd.read_csv('tmp.tbl', sep='\s+', skiprows=3, comment='#', usecols=[0, 2], header=None,
	                          names=['protein_id', 'pVOG'])
	os.remove('tmp.tbl')
	return hmm_results


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
		merged_df = pd.read_csv(outfile, sep='\t')
	else:
		logger.info('Calling genes with prodigal')
		cmd = f"""prodigal -a tmp.fa -f gff -i {infile} -m -p meta -q"""
		stdout, stderr = execute(cmd)
		if stderr:
			logger.info(stderr.decode('utf-8'))
		if stdout:
			logger.debug(stdout.decode('utf-8'))

		rgx = re.compile('((\S+)_(\S+))\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+([0-9\-]+).*')
		protein_results = []
		cleaned_proteins = []
		with open('tmp.fa', 'r') as handle:
			for record in SeqIO.parse(handle, 'fasta'):
				m = rgx.search(record.description)
				if m:
					protein_results.append((m.group(2), m.group(1), m.group(4), m.group(5), m.group(6), 'prodigal'))
					record.id = f'{m.group(2)}_{m.group(1)}'
					record.name = ''
					record.description = ''
					record.seq = record.seq[:-1]
					cleaned_proteins.append(record)
		SeqIO.write(cleaned_proteins, f'{args.output_folder}/prodigal.fa', 'fasta')
		hmm_results = runHMMScan(f'{args.output_folder}/prodigal.fa')
		hmm_results['contig'], hmm_results['gene_id'] = hmm_results['protein_id'].str.split('_', 1).str
		hmm_results = hmm_results[['contig', 'gene_id', 'pVOG']]
		df = pd.DataFrame(protein_results, columns=['contig', 'gene_id', 'start', 'end', 'frame', 'method'])
		merged_df = df.merge(hmm_results, how='left', on=['contig', 'gene_id'])
		logging.info(f'Writing file {os.path.abspath(outfile)}')
		merged_df.to_csv(os.path.abspath(outfile), sep='\t', index=False)
		os.remove('tmp.fa')

	return merged_df


if __name__ == "__main__":
	main()
