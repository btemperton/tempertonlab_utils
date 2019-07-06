import argparse
import glob
import logging
import os
import subprocess
from Bio import SeqIO
from string import Template
import gzip

template = """
#PBS -V # set verbose output
#PBS -d . # set working directory to .
#PBS -q pq # submit to the pq (max 96 cores)
#PBS -n map.$SAMPLE
#PBS -l walltime=12:00:00
#PBS -A Research_Project-172179 # research project to submit under.
#PBS -l nodes=1:ppn=16 # or nodes=number of nodes required. ppn=number of processors per node
#PBS -j oe

cd "/local/pbstmp.$${PBS_JOBID}"

FVE_LOCATION=$FVE_LOCATION

source activate $CONDA_ENV

COUNT=$$(zcat $FWD_READS | awk '{s++}END{print s/4}')

if [ "$$COUNT" -ge $SUBSAMPLE ]; then

	seqtk sample -s 42 $FWD_READS $SUBSAMPLE > fwd.fq
	seqtk sample -s 42 $REV_READS $SUBSAMPLE > rev.fq

	java -Xmx108G -cp $$FVE_LOCATION/bin FastViromeExplorer \
	-1 fwd.fq \
	-2 rev.fq \
	-i $KALLISTO_INDEX \
	-l $VIRAL_LIST \
	-cr $MIN_RATIO \
	-co $MIN_COVERAGE \
	-cn $MIN_READS \
	-reportRatio true


	cp log.txt $OUTPUT_FOLDER/$SAMPLE.fve.log
	cp FastViromeExplorer-final-sorted-abundance.tsv $OUTPUT_FOLDER/$SAMPLE.fve.sorted.abundance.tsv
	cp abundance.tsv $OUTPUT_FOLDER/$SAMPLE.fve.abundance.tsv
fi


"""

def main():
	parser = argparse.ArgumentParser(description="""""")
	parser.add_argument('--metagenomes', '-m', dest='metagenomes', required=True, help='Comma-separated file of metagenomes in the format SAMPLE_NAME,Forward reads, Reverse Reads (no header)')
	parser.add_argument('--fasta', '-f', dest='fasta_file', required=True, help='Fasta file of the reference sequences')
	parser.add_argument('--output_dir', '-d', dest='output_folder', required=True, help='Output directory')
	parser.add_argument('--subsample_no', '-S', dest='subsample' ,type=int, default=10000000)
	parser.add_argument('--log', '-l', dest='logfile', default='create_kallisto.index.log')
	parser.add_argument('--fve_loc', dest='fve_loc', required=True, help='The location of the FastViromeExplorer directory')
	parser.add_argument('--min_ratio', dest='min_ratio', type=float, default=0.3, help='Minimum ratio (see https://peerj.com/articles/4227/)')
	parser.add_argument('--min_coverage', dest='min_coverage', type=float, default=0.1, help='Minimum coverage (see https://peerj.com/articles/4227/)')
	parser.add_argument('--min_reads', dest='min_reads', type=int, default=10, help='Minimum number of reads (see https://peerj.com/articles/4227/)')
	parser.add_argument('--overwrite', dest='overwrite', type=bool, default=False)
	parser.add_argument('--conda_env', dest='conda_env', default='fve')

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
	create_kallisto_index()


	with open(args.metagenomes, 'r') as handle:
		for line in handle.readlines():
			bits = line.strip().split(',')
			launch_job(bits[0], bits[1], bits[2])

def launch_job(sample, fwd, rev):

	if not os.path.isfile('{}/{}.fve.sorted.abundance.tsv'.format(args.output_folder, sample)) or args.overwrite:
		d = {'SAMPLE': sample,
			 'OUTPUT_FOLDER': os.path.abspath(args.output_folder),
			 'FWD_READS': os.path.abspath(fwd),
			 'REV_READS': os.path.abspath(rev),
			 'FVE_LOCATION': os.path.abspath(args.fve_loc),
			 'MIN_RATIO' : args.min_ratio,
			 'MIN_COVERAGE': args.min_coverage,
			 'MIN_READS': args.min_reads,
			 'SUBSAMPLE': args.subsample,
			 'KALLISTO_INDEX': f'{os.path.abspath(args.output_folder)}/kallisto-index.idx',
			 'VIRAL_LIST': f'{os.path.abspath(args.output_folder)}/fve_index_file.txt',
			 'CONDA_ENV': args.conda_env}
		t = Template(template).substitute(d)
		logging.debug('Writing job file to {}/{}.job'.format(args.output_folder, sample))

		with open('{}/{}.job'.format(args.output_folder, sample), 'w') as handle:
			handle.write(t)
		cmd = 'qsub {}/{}.job'.format(args.output_folder, sample)
		stdout, stderr = execute(cmd)
		logger.info(stderr)
		logger.debug(stdout)
	else:
		logger.debug('No need to launch job for {} as output file already exists'.format(sample))

def create_kallisto_index():
	if not os.path.isfile(f'{os.path.abspath(args.output_folder)}/fve_index_file.txt'):
		cmd = f"""bash {args.fve_loc}/utility-scripts/generateGenomeList.sh {os.path.abspath(args.fasta_file)} {os.path.abspath(args.output_folder)}/fve_index_file.txt"""
		execute(cmd)
	else:
		logging.info(f'Index file containing contig lengths already exists at {os.path.abspath(args.output_folder)}/fve_index_file.txt')

	if not os.path.isfile(f'{os.path.abspath(args.output_folder)}/kallisto-index.idx'):
		cmd = f'kallisto index -i {os.path.abspath(args.output_folder)}/kallisto-index.idx {os.path.abspath(args.fasta_file)}'
		execute(cmd)
	else:
		logging.info(f'Kallisto index file already exists at {os.path.abspath(args.output_folder)}/kallisto-index.idx')


def create_output_dir():
	try:
		os.mkdir(args.output_folder)
		logger.debug('Created output folder at {}'.format(args.output_folder))
	except OSError:
		logger.debug('Output folder at {} already exists'.format(args.output_folder))

def execute(command):
	logger.info('\nExecuting {}'.format(command))
	process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
							   shell=True)
	(stdout, stderr) = process.communicate()

	return stdout, stderr


if __name__ == "__main__":
	main()
