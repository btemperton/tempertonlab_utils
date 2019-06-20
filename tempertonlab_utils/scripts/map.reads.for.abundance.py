import argparse
import glob
import logging
import os
import subprocess
from string import Template

template = """
#PBS -V # set verbose output
#PBS -d . # set working directory to .
#PBS -q pq # submit to the pq (max 96 cores)
#PBS -n map.$SAMPLE
#PBS -l walltime=08:00:00
#PBS -A Research_Project-172179 # research project to submit under.
#PBS -l nodes=1:ppn=$THREADS # or nodes=number of nodes required. ppn=number of processors per node
#PBS -j oe

cd "/local/pbstmp.$${PBS_JOBID}"

##This assumes a conda environment has been created with the following command:
## conda create -y -n calculate.viral.abundance bwa samtools bowtie2 bamm pysam bbmap

source activate calculate.viral.abundance

BASECAMP=$OUTPUT_FOLDER
SAMPLE_NAME=$SAMPLE
BOWTIE_INDEX=$$BASECAMP/ref
FWD_READS=$FWD_READS
REV_READS=$REV_READS
LOGFILES=$$BASECAMP/logs
THREADS=$THREADS
OUTPUT_DIR=$$BASECAMP

mkdir -p $$LOGFILES

if [[ $$FWD_READS == ftp* ]]
then
    wget -O tmp.fwd.fq.gz $$FWD_READS
else
    cp $$FWD_READS tmp.fwd.fq.gz
fi

if [[ $$REV_READS == ftp* ]]
then
    wget -O tmp.rev.fq.gz $$REV_READS
else
    cp $$REV_READS tmp.rev.fq.gz
fi

clumpify.sh -da in1=tmp.fwd.fq.gz in2=tmp.rev.fq.gz out=clumped.fq.gz dedupe optical
ln -s clumped.fq.gz temp.fq.gz

filterbytile.sh -da in=temp.fq.gz out=filtered_by_tile.fq.gz
rm temp.fq.gz; ln -s filtered_by_tile.fq.gz temp.fq.gz

bbduk.sh -da in=temp.fq.gz out=trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=50 ref=adapters ftm=5 ordered
rm temp.fq.gz; ln -s trimmed.fq.gz temp.fq.gz

bbduk.sh -da in=temp.fq.gz out=filtered.fq.gz k=31 ref=artifacts,phix ordered cardinality
rm temp.fq.gz; ln -s filtered.fq.gz temp.fq.gz

#converts it back into paired files (not interleaved) and subsamples down to 10m reads
reformat.sh -da sampleseed=42 samplereadstarget=10000000 in=temp.fq.gz out1=fwd.fq.gz out2=rev.fq.gz

bowtie2 -x $$BOWTIE_INDEX \
-1 fwd.fq.gz \
-2 rev.fq.gz \
-S mapping.sam \
--threads $$THREADS 2>&1 | tee $$LOGFILES/$$SAMPLE_NAME.bowtie2.mapping.log

samtools view --threads $$THREADS -F 4 -bS -o mapping.bam mapping.sam
samtools sort --threads $$THREADS -o mapping.sorted.bam mapping.bam
samtools index -@ $$THREADS mapping.sorted.bam

bamm filter -b mapping.sorted.bam --percentage_id 0.95 --percentage_aln 0.9
mv mapping.sorted_filtered.bam $$OUTPUT_DIR/$$SAMPLE_NAME.bamm.id95.aln90.bam
mv mapping.sorted_filtered.bam.bai $$OUTPUT_DIR/$$SAMPLE_NAME.bamm.id95.aln90.bam.bai

"""


def main():
	parser = argparse.ArgumentParser(description='Create mapping files for viral genomes against a set of metagenomes')
	parser.add_argument('--fasta', '-f', dest='fasta_file', required=True, help='Fasta file of the reference sequences')
	parser.add_argument('--metagenomes', '-m',
	                    dest='metagenomes',
	                    required=True,
	                    help='Comma-separated file of metagenomes in the format SAMPLE_NAME,Forward reads, Reverse Reads (no header)')
	parser.add_argument('--output', '-o', dest='output_folder', required=True,
	                    help='Output folder where the files will be created (will be generated if it doesn\'t exist)')
	parser.add_argument('--threads', '-t', dest='threads', default=16, type=int,
	                    help='Number of threads to use')
	parser.add_argument('--log', '-l', dest='logfile', default='map.reads.for.abundance.log')
	parser.add_argument('--overwrite', dest='overwrite', type=bool, default=False)
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
	create_bowtie2_index()

	with open(args.metagenomes, 'r') as handle:
		for line in handle.readlines():
			bits = line.strip().split(',')
			launch_job(bits[0], bits[1], bits[2])


def launch_job(sample, fwd, rev):
	if not os.path.isfile('{}/{}.bamm.id95.aln90.bam'.format(args.output_folder, sample)) or args.overwrite:
		d = {'SAMPLE': sample,
		     'THREADS': args.threads,
		     'OUTPUT_FOLDER': os.path.abspath(args.output_folder),
		     'FWD_READS': fwd,
		     'REV_READS': rev}
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


def create_output_dir():
	try:
		os.mkdir(args.output_folder)
		logger.debug('Created output folder at {}'.format(args.output_folder))
	except OSError:
		logger.debug('Output folder at {} already exists'.format(args.output_folder))


def create_bowtie2_index():
	if len(glob.glob('{}/ref*.bt2'.format(args.output_folder))) != 6:
		cmd = "bowtie2-build {} {}/ref".format(args.fasta_file, args.output_folder)
		stdout, stderr = execute(cmd)
		logger.info(stderr)
		logger.debug(stdout)
	else:
		logger.info('Bowtie2 index already exists at {}/ref'.format(args.output_folder))


def execute(command):
	logger.info('\nExecuting {}'.format(command))
	process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
	                           shell=True)
	(stdout, stderr) = process.communicate()

	return stdout, stderr


if __name__ == "__main__":
	main()
