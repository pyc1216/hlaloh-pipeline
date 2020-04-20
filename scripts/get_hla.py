#!/usr/bin/python
# -*- coding: utf-8 -*-  


##################################################################################
###input：OptiType result tsv [PolySolver result counts1.R0k6 & counts2.R0k6     #
###output: lohhla input hlatype file (--hlaPath )                                #
###purpose:                                                                      #
###                                                                              #
###                                                                              #
##################################################################################
                                                                                  
################################ input format ####################################
#$ cat normal_result.tsv                                                         #
#	A1	A2	B1	B2	C1	C2	Reads	Objective                                    #
#0	A*02:01	A*11:01	B*27:05	B*54:01	C*01:02	C*01:02	13374	12892.535999999987   #
##################################################################################
                                                                                  
################################ output format ###################################
#$ cat loh/hla/winners.hla.txt                                                   #
#HLA-A	hla_a_11_01_01	hla_a_11_01_01                                           #
#HLA-B	hla_b_54_01_03	hla_b_27_05_02                                           #
#HLA-C	hla_c_01_02_01	hla_c_01_02_01                                           #
##################################################################################

from __future__ import print_function
import sys
import os


def csv_parser(input_opti):
	list_hla = []
	fmt_converter = lambda x: x.lower().replace('-','_').replace('*', '_').replace(':', '_')
	with open(input_opti) as handle:
		lines = handle.readlines()
	if len(lines) == 2:
		list_hla.extend(['hla_' + fmt_converter(x) for x in lines[1].split('\t')[1:7]])
	return list_hla

def write_output(list_hla, output):
	hla_genes = ['hla_a', 'hla_b', 'hla_c']
	with open(output, 'w') as handle:
		if not list_hla: #empty list_hla
			print("[WARNING] Nothing was written to {}!".format(output))
			return None
		for gene in hla_genes:
			try: #if missing gene,then continue
				allele1, allele2 = [x for x in list_hla if x.startswith(gene)]
			except ValueError:
				continue
			if allele1 == allele2: #only one uniq allele
				print('[INFO] {} is HOM, allele={}'.format(gene, allele1))
			else:
				alleles = allele1 + '\t' + allele2
				line = gene.replace('_', '-').upper() + '\t' + alleles + '\n'
				handle.write(line)

def counts_parser(counts):
	zip_allele_scores = []
	with open(counts) as handle:
		for line in handle:
			allele, score = line.strip().split('\t')
			score = float(score)
			zip_allele_scores.append((allele, score))
	return sorted(zip_allele_scores, key=lambda x: x[1], reverse=True)

def get_rank_score(zip_allele_scores_gene, allele):
	#zip_allele_scores_gene = [x for x in zip_allele_scores if x[0].startswith(gene)]
	rank = len(zip_allele_scores_gene) + 1
	score = None
	allele_poly = allele
	for i, item in enumerate(zip_allele_scores_gene):
		allele_poly, score_poly = item
		score_poly = float(score_poly)
		if i == 0:
			max_score = score_poly
		if allele_poly.startswith(allele):
			rank = i + 1
			score_percent = score_poly / max_score
			break
	return (allele_poly, rank, score_percent)

def redefine_poly_allele1(zip_allele_scores_gene, allele1, allele2):
	#zip_allele_scores_gene = [x for x in zip_allele_scores if x[0].startswith(gene)]
	num_alleles = len(zip_allele_scores_gene) + 1
	metrix_1 = [allele1, num_alleles, None] #[allele1, rank1, score1]
	metrix_2 = [allele2, num_alleles, None]
	for i, item in enumerate(zip_allele_scores_gene):
		allele_poly, score_poly = item #str, float
		if i == 0:
			max_score = score_poly
		if metrix_1[1] != num_alleles and metrix_2[1] != num_alleles:
			break
		if allele_poly.startswith(allele1):
			#rank1 = i + 1
			#score_percent1 = score_poly / max_score
			metrix_1 = [allele_poly, i + 1, score_poly / max_score]
		if allele_poly.startswith(allele2):
			metrix_2 = [allele_poly, i + 1, score_poly / max_score]
	metrix_allele = metrix_1 if metrix_1[1] < metrix_2[1] else metrix_2
	#High confidence allele1: rank < 50 and percent of best max score > 0.8
	if metrix_allele[1] < 50 and metrix_allele[2] > 0.8:
		return metrix_allele[0]
	#Low confidence allele1:
	elif metrix_allele[1] < 100 and metrix_allele[2] > 0.5:
		print("[WARNING] Low confidence allele1! Still accept. {}(rank={:.2f}, score_percent={:.2f})".format(metrix_allele[0], metrix_allele[1], metrix_allele[2]))
		return metrix_allele[0]
	#Uncertain situation:
	else:
		print("[WARNING] Uncertain allele1! Drop. {}(rank={:.2f}, score_percent={:.2f}).".format(metrix_allele[0], metrix_allele[1], metrix_allele[2]))
		return None
	#return metrix_allele[0] if metrix_allele[1] < 50 and metrix_allele[2] > 0.8 else None


def main():
	if len(args) == 2:
		input_opti, output = sys.argv[1:]
		list_hla = csv_parser(input_opti) #HLA-A*11:01,HLA-A*02:01,HLA-B*27:05,HLA-B*54:01,HLA-C*01:02
		write_output(list_hla, output)
	else:
		input_opti, input_poly, output = args
		list_hla = csv_parser(input_opti)
		counts1 = os.path.join(input_poly, 'counts1.R0k6')
		counts2 = os.path.join(input_poly, 'counts2.R0k6')
		if not os.path.exists(counts1) or not os.path.exists(counts2):
			print("[WARNING] Missing PolySolver counts file: {} and {}".format(counts1, counts2))
			write_output(list_hla, output)
		else:
			list_hla_out = []
			zip_allele1_scores = counts_parser(counts1)
			zip_allele2_scores = counts_parser(counts2)
			hla_genes = ['hla_a', 'hla_b', 'hla_c']
			for gene in hla_genes:
				zip_allele1_scores_gene = [x for x in zip_allele1_scores if x[0].startswith(gene)]
				zip_allele2_scores_gene = [x for x in zip_allele2_scores if x[0].startswith(gene)]
				allele1, allele2 = [x for x in list_hla if gene in x]
				poly_allele1 = zip_allele1_scores_gene[0][0]
				
				if allele1 == allele2:
					print("[INFO] {} is HOM, allele={}, skipping this allele".format(gene, allele1))
					continue
				
				#if OptiType's allele1&2 dose not match PolySolver's allele1 either, so re-define poly_allele1
				if not poly_allele1.startswith(allele1) and not poly_allele1.startswith(allele2):
					print("[INFO] OptiType two alleles({}, {}) do not matches PolySolver allele1({}). Try to get another allele1 in {}..."\
					.format(allele1, allele2, poly_allele1, counts1))
					poly_allele1 = redefine_poly_allele1(zip_allele1_scores_gene, allele1, allele2)
					if not poly_allele1:
						#print("OptiType two alleles({}, {}) not matches PolySolver allele1({}) !".format(allele1, allele2, poly_allele1))
						print("[WARNING] Failed in trying to get anthoer allele1.")
						continue
				
				if poly_allele1.startswith(allele1):
					pass
				elif poly_allele1.startswith(allele2):
					allele1, allele2 = allele2, allele1 #swap allele1 and allele2:

				poly_allele2, rank, score_percent = get_rank_score(zip_allele2_scores_gene, allele2)
				#High confidence allele2
				if rank < 100 and score_percent > 0.5:
					list_hla_out.extend([poly_allele1, poly_allele2])
				#Low confidance allele2
				else:
					print("[WARNING] OptiType allele2({}) not be supported by PolySolver ({}, rank={:.2f}, score_percent={:.2f}) !".format(allele2, poly_allele2, rank, score_percent))
					list_hla_out.extend([poly_allele1, poly_allele2])

			write_output(list_hla_out, output)


if __name__ == '__main__':
	args = sys.argv[1:]
	if len(args) < 2 or len(args) > 3:
		print('[USAGE] python {} normal_result.tsv(OptiType) [ loh/polysolver/(PolySolver dir) ] output(hla_optitype) '.format(sys.argv[0]))
		exit(1)
	else:
		main()
