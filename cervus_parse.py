#!/usr/bin/env python

import sys, zipfile, os, shutil, vcf, argparse, re
from collections import defaultdict
from collections import OrderedDict
from gooey import Gooey, GooeyParser

class Parser():
	def __init__(self, tempdir, infile, outfile):
		self.infile = infile
		self.outfile = outfile if os.path.splitext(outfile)[1] == '.csv' else (outfile + '.csv')
		self.flags = defaultdict(dict)
		vcf_files = self.extract_data(tempdir)
		self.hotspot_d = self.get_hotspot_info(vcf_files[0])
		self.data_d = self.read_vcfs(vcf_files)
		self.sample_l = self.get_sample_order(vcf_files)

	
	def extract_data(self, tempdir):
		if os.path.isdir(tempdir) == True:
			shutil.rmtree(tempdir)
		os.mkdir(tempdir)
		with zipfile.ZipFile(self.infile) as z:
			z.extractall(tempdir)
		return [ os.path.join(tempdir, file) for file in os.listdir(tempdir) if file.endswith('.vcf.gz') ]


	def assign_genotype(self, vcf_record, sample_name):
		genotype = []
		flagged = False
		for allele in vcf_record.samples[0]['GT'].split('/'):
			if allele == '.' or len(vcf_record.FILTER) > 0:
				self.report_flag(vcf_record, sample_name, 'No Call')
				# This is where you need to flag No Call genotype is present
				return [ '0', '0' ]
			elif int(allele) > 1:
				genotype.append('2')
				if flagged == False:
					self.report_flag(vcf_record, sample_name, 'Unusual genotype: {}'.format(vcf_record.samples[0]['GT']))
					flagged = True
				# This is where you need to flag when a unusual genotype is present
			else:
				genotype.append(str(int(allele)+1))
		return genotype


	def report_flag(self, vcf_record, sample_name, flag):
		chrom = vcf_record.CHROM
		pos = vcf_record.POS
		id = vcf_record.ID
		ref = vcf_record.REF
		alts = vcf_record.ALT
		flag = flag
		filters = vcf_record.FILTER
		qual = vcf_record.QUAL
		types = vcf_record.INFO['TYPE']
		if (chrom, str(pos)) in self.flags[sample_name]:
			self.flags[sample_name][(chrom, str(pos))][5].append(flag)
		else:
			self.flags[sample_name][(chrom, str(pos))] = [chrom, str(pos), id, ref, alts, [flag], filters, str(qual), types]
	
	def check_types_return_true_if_wrong_base(self, vcf_record, sample_name):
		called_snps = 0
		for type in vcf_record.INFO['TYPE']:
			if type == 'snp':
				called_snps += 1
			else:
				self.report_flag(vcf_record, sample_name, 'Complex/MNP')
		if called_snps > 1:
			self.report_flag(vcf_record, sample_name, 'Wrong hotspot base')
			return True
		return False
	
	def check_borderline_quality(self, vcf_record, sample_name):
		if vcf_record.QUAL >= 10 and vcf_record.QUAL < 100:
			self.report_flag(vcf_record, sample_name, 'Borderline quality')
	
	# def check_correct_alt_base(self, vcf_record, sample_name):
		# for alt in vcf_record.ALT:
			# if not str(alt) in vcf_record.ALT
	
	
	def make_flags_all_str(self):
		for sample_name in self.flags.keys():
			for chrom, pos in self.flags[sample_name].keys():
				[chrom1, pos1, id1, ref1, alts1, flags1, filters1, qual1, types1] = self.flags[sample_name][(chrom, str(pos))]
				alts2 = '|'.join(str(x) for x in alts1)
				flags2 = '|'.join(str(x) for x in flags1)
				filters2 = '|'.join(str(x) for x in filters1)
				types2 = '|'.join(str(x) for x in types1)
				self.flags[sample_name][(chrom, str(pos))] = [chrom1, pos1, id1, ref1, alts2, flags2, filters2, qual1, types2]
	
	
	def read_vcfs(self, vcf_files):
		data_d = defaultdict(dict)
		for file in vcf_files:
			vcf_reader = vcf.Reader(filename=file)
			sample_name = ''
			for rec in vcf_reader:
				if sample_name == '':
					sample_name = rec.samples[0].sample
				# This would be the place to add flags for wrong mutation (compare ALT and OALT?) and flag when hotspot is part of a MNP
				try:
					if rec.INFO['HS'] == True:
						self.check_borderline_quality(rec, sample_name)
						if self.check_types_return_true_if_wrong_base(rec, sample_name) == True:
							data_d[sample_name][rec.ID] = [ '0', '0' ]
						else:
							data_d[sample_name][rec.ID] = self.assign_genotype(rec, sample_name)
					else:
						print("Error: HS flag = {}".format(rec.INFO['HS']))
						sys.exit(1)
				except KeyError:
					continue
		return data_d


	def get_hotspot_info(self, vcf_file):
		hotspot_d = OrderedDict()
		vcf_reader = vcf.Reader(filename=vcf_file)
		for rec in vcf_reader:
			try:
				if rec.INFO['HS'] == True:
					for i, oid in enumerate(rec.INFO['OID']):
						if oid != None and not (oid in hotspot_d):
							hotspot_d[oid] = (rec.INFO['OPOS'][i], rec.INFO['OREF'][i], rec.INFO['OALT'][i])
				else:
					print("Error: HS flag = {}".format(rec.INFO['HS']))
					sys.exit(1)
			except KeyError:
				continue
		return hotspot_d


	def get_sample_order(self, vcf_files):
		sample_l = []
		for file in vcf_files:
			vcf_reader = vcf.Reader(filename=file)
			sample_name = ''
			for rec in vcf_reader:
				if sample_name == '':
					sample_l.append(rec.samples[0].sample)
					break
		return sample_l


	def output_tsv(self, data_d, hotspot_d, sample_l):
		out_l = ['ID']
		outfile = open(self.outfile, 'w')
		outflags = open("{}.flags.csv".format(self.outfile), 'w')
		for k, v in hotspot_d.items():
			opos, oref, oalt = v
			out_l.append(str(opos))
			out_l.append(str(opos))
		outfile.write('{}\n'.format(','.join(out_l)))
		self.make_flags_all_str()
		for sample in sample_l:
			out_l = [sample]
			print("Writing sample {}".format(sample))
			for k, v in hotspot_d.items():
				out_l += data_d[sample][k]
			outfile.write('{}\n'.format(','.join(out_l)))
			outflags.write('Sample: {}\n'.format(sample))
			outflags.write('Chrom,Position,ID,Ref,Variant,Flag,Filter,Quality,Type\n')
			for flag in sort_nicely(self.flags[sample].values()):
				outflags.write('{}\n'.format(','.join(flag)))
			outflags.write('\n')
		outfile.close()


def sort_nicely(l):
	""" Sort the given list in the way that humans expect.
	"""
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key[0]) ]
	return sorted(l, key=alphanum_key )


@Gooey(default_size=(1000, 530))
def main():
	tempdir = '.tmp'
	parser = GooeyParser(description="Cervus Parser - Convert vcf.zip files from Ion Torrent to Cervus input.")
	parser.add_argument("Infile", help="VCF.ZIP file downloaded from Ion Torrent.", action="store", widget='FileChooser')
	parser.add_argument("Outfile", help="Desired name for your Cervus-formatted file. This will be added to the same folder as the input file.", action="store")
	args = parser.parse_args()
	cervus_parser = Parser(tempdir, args.Infile, os.path.join(os.path.dirname(args.Infile), args.Outfile))
	cervus_parser.output_tsv(cervus_parser.data_d, cervus_parser.hotspot_d, cervus_parser.sample_l)
	shutil.rmtree(tempdir)


if __name__ == '__main__':
	main()
