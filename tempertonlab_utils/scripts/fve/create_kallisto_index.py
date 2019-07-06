import argparse
import logging
import os
from string import Template
import gzip

template = """
#PBS -V # set verbose output
#PBS -d . # set working directory to .
#PBS -q highmem # submit to the pq (max 96 cores)
#PBS -n make.kallisto
#PBS -l walltime=12:00:00
#PBS -A Research_Project-172179 # research project to submit under.
#PBS -l nodes=1:ppn=32 # or nodes=number of nodes required. ppn=number of processors per node
#PBS -j oe

cd "/local/pbstmp.$${PBS_JOBID}"

FVE_LOCATION=$FVE_LOCATION

source activate $CONDA_ENV
cp $INPUT_FILE tmp.fa.gz
unpigz tmp.fa.gz

bash $$FVE_LOCATION/utility-scripts/generateGenomeList.sh tmp.fa fve_index_file.txt

kallisto index -i tmp.fa kallisto-index.idx


cp fve_index_file.txt $OUTPUT_DIR
cp kallisto-index.idx $OUTPUT_DIR
"""

def main():
	parser = argparse.ArgumentParser(description="""Create Kallisto index for viral genomes against a set of metagenomes using FastViromeExplorer. It is sensible to do this if you are running lots of samples against the same database.

    You will need to have set up a conda environment for running this as follows:

    conda create -n fve samtools kallisto=0.43.1 openjdk python=3.7 ipython biopython pigz numpy pandas seqtk

    This will also be used for the running of FastViromeExplorer.
    """)
	parser.add_argument('--output_dir', '-d', dest='output_dir', required=True, help='Location for the kallisto index and files')
    parser.add_argument('--input', , '-i', dest='input_file', required=True, help='Gzipped fasta file to be indexed for FVE')
    parser.add_argument('--prefix', , '-p', dest='prefix', required=True, help='Prefix for the output file names')
	parser.add_argument('--log', '-l', dest='logfile', default='create_kallisto.index.log')
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



def create_output_dir():
	try:
		os.mkdir(args.output_folder)
		logger.debug('Created output folder at {}'.format(args.output_folder))
	except OSError:
		logger.debug('Output folder at {} already exists'.format(args.output_folder))



if __name__ == "__main__":
	main()
