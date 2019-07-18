import argparse
import gzip
import logging

import pandas as pd
import pysam
from Bio import SeqIO


def main():
	parser = argparse.ArgumentParser(description='Identify HVRs in a genome')
	parser.add_argument('--bam', '-b', dest='bam_file', required=True, help='bam file from the mapping')
	parser.add_argument('--output', '-o', dest='output_file', required=True,
	                    help='Output file listing the reads mapped')
	parser.add_argument('--fwd', '-1', dest='fwd_reads', required=True)
	parser.add_argument('--rev', '-2', dest='rev_reads', required=True)
	parser.add_argument('--sample', '-s', dest='sample_name', required=True)
	parser.add_argument('--cutoff', '-c', dest='cutoff', required=True)
	parser.add_argument('--log', '-l', dest='logfile', default='parse.bam.for.RPKG.log')

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



	read_count, gbp = count_reads()
	df = parse_bam_for_RPKG()
	df['sample'] = args.sample_name
	df['cutoff'] = args.cutoff
	df['read_count_million'] = read_count
	df['sample_gbp'] = gbp
	df['RPKG'] = df['read_count'] / df['length_kbp'] / gbp
	df.to_csv(args.output, sep='\t', index=False)

def count_reads():
	read_count = 0
	read_length = 0
	with gzip.open(args.fwd_reads, 'rt') as handle:
		for record in SeqIO.parse(handle, 'fastq'):
			read_count +=1
			read_length += len(record.seq)
	with gzip.open(args.rev_reads, 'rt') as handle:
		for record in SeqIO.parse(handle, 'fastq'):
			read_count +=1
			read_length += len(record.seq)

	return read_count/1e6, read_length/1e9


def parse_bam_for_RPKG():
	samfile = pysam.AlignmentFile(args.bam_file, "rb")
	refs = [(x, samfile.get_reference_length(x)/1000,
	         samfile.count(contig=x, read_callback='all')) for x in samfile.references]
	return pd.DataFrame(refs, columns = ['contig', 'length_kbp', 'reads_mapped'])


if __name__ == "__main__":
	main()