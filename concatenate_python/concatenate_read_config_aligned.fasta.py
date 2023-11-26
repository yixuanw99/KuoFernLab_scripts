#! /home/lykuo/miniconda2/bin python

import os
import collections
import sys

#Define parameter name
import_align_dir=sys.argv[1]
TaxaName=sys.argv[2]

coding_list_length=int(sys.argv[3])
#noncoding_list_length=int(sys.argv[4])

coding_list_end=5+coding_list_length
#noncoding_list_end=coding_list_end+noncoding_list_length

coding_list=sys.argv[5:coding_list_end]
#noncoding_list=sys.argv[coding_list_end:noncoding_list_end]

alignment_FileName = TaxaName + "_cpDNA" + ".fasta"
presence_table_FileName = TaxaName + "_cpDNA" + ".csv"
parition_file_FileName = TaxaName + "_cpDNA_partition" + ".txt"
parition12_file_FileName = TaxaName + "_cpDNA_partition_codon12" + ".txt"
parition3_file_FileName = TaxaName + "_cpDNA_partition_codon3" + ".txt"

#--------------------
# import_align_dir = '/home2/analyses/Lomariopsidaceae/cpDNA_alignment/20221004/'
output_alignment = open(alignment_FileName, "w")
output_presence_table = open(presence_table_FileName, "w")
output_parition_file = open(parition_file_FileName, "w")
output_parition12_file = open(parition12_file_FileName, "w")
output_parition3_file = open(parition3_file_FileName, "w")


taxa_to_fastaID_seq_dic = {}
loci_taxa_presence_dict = {}
cumulative_len = 0
fasta_list = []
taxa_raw_list = []

#phylip to fasta#
#phylip_list = ["L100", "L101", "L102", "L103", "L106", "L107", "L108", "L112", "L113", "L114", "L115", "L116", "L117", "L118", "L119", "L120", "L121", "L122", "L123", "L124", "L125", "L126", "L127", "L129", "L130", "L131", "L132", "L133", "L135", "L137", "L138", "L139", "L141", "L142", "L143", "L144", "L145", "L146", "L147", "L148", "L149", "L150", "L151", "L152", "L153", "L154", "L155", "L156", "L157", "L158", "L159", "L160", "L161", "L162", "L163", "L164", "L166", "L167", "L168", "L169", "L170", "L171", "L172", "L173", "L174", "L175", "L176", "L177", "L178", "L179", "L180", "L181", "L182", "L183", "L184", "L185", "L186", "L188", "L189", "L190", "L191", "L193", "L194", "L196", "L197", "L198", "L199", "L200", "L201", "L202", "L203", "L204", "L207", "L209", "L210", "L211", "L212", "L213", "L214", "L215", "L217", "L218", "L219", "L220", "L221", "L222", "L223", "L224", "L225", "L226", "L227", "L228", "L229", "L230", "L231", "L232", "L234", "L235", "L236", "L237", "L238", "L239", "L240", "L241", "L242", "L243", "L245", "L246", "L247", "L248", "L249", "L250", "L252", "L253", "L254", "L255", "L256", "L257", "L258", "L260", "L261", "L262", "L263", "L264", "L265", "L266", "L267", "L268", "L269", "L270", "L271", "L272", "L273", "L274", "L275", "L276", "L277", "L278", "L279", "L280", "L281", "L282", "L283", "L284", "L285", "L286", "L289", "L290", "L291", "L292", "L293", "L294", "L295", "L296", "L297", "L298", "L299", "L300", "L302", "L303", "L304", "L305", "L306", "L307", "L308", "L309", "L310", "L311", "L312", "L313", "L314", "L315", "L316", "L317", "L318", "L319", "L320", "L321", "L322", "L323", "L324", "L325", "L326", "L327", "L328", "L329", "L331", "L332", "L333", "L334", "L335", "L336", "L337", "L338", "L339", "L340", "L341", "L342", "L343", "L344", "L345", "L346", "L347", "L348", "L349", "L350", "L351", "L352", "L353", "L354", "L355", "L356", "L357", "L358", "L359", "L365", "L366", "L367", "L368", "L369", "L370", "L372", "L373", "L374", "L375", "L376", "L377", "L379", "L380", "L381", "L382", "L383", "L384", "L385", "L386", "L387", "L388", "L389", "L391", "L392", "L394", "L395", "L396", "L397", "L398", "L400", "L401", "L402", "L403", "L404", "L405", "L406", "L407", "L408", "L409", "L410", "L411", "L412", "L413", "L415", "L416", "L417", "L418", "L419", "L420", "L421", "L422", "L423", "L424", "L425", "L426", "L428", "L429", "L430", "L431", "L432", "L433", "L434", "L435", "L436", "L437", "L438", "L439", "L440", "L441", "L442", "L443", "L447", "L448", "L449", "L450", "L85", "L90", "L91", "L92", "L93", "L94", "L95", "L97", "L98", "L99"]
#phylip_dir = '/home/lykuo/lab_data/NGS_data/GoFlag/Marattiales/Marattiales/Marattiales_Target/'
#for phylip in phylip_list:
#	phylip_file = phylip_dir + phylip + '.keep1.prune.phy'
#	output_fasta = open(phylip_dir + phylip + '.keep1_prune_fasta', "w")
#	with open(phylip_file, "r") as f:
#		for row in f.readlines()[1:]:
#			header = row.split("\t")[0]
#			seq = row.split("\t")[1]
#			fas = '>' + header + '\n' + seq
#			output_fasta.write(fas)
#	output_fasta.close

#read fasta#
for file in os.listdir(import_align_dir):
	name_length = len(file.split('.'))
	if file.split('.')[name_length-1] == 'fas': #file type filelsit and taxalist
		fasta_list.append(file)
		fasta_file_name = import_align_dir + file
		for line in open(fasta_file_name, "r"):
			if line.startswith('>'):
				line = line.rstrip()
				fa_ID = line.replace('>','')
				taxa_raw_list.append(fa_ID)

taxa_list = list(set(taxa_raw_list))
#fasta_list = sorted(fasta_list)

for taxon in taxa_list:
	taxa_to_fastaID_seq_dic[taxon] = ""

#locus_list_cp = ["accD", "atpA", "atpB", "atpE", "atpF", "atpH", "atpI", "ccsA", "cemA", "chlB", "chlL", "chlN", "clpP", "infA", "matK", "ndhA", "ndhB", "ndhC", "ndhD", "ndhE", "ndhF", "ndhG", "ndhH", "ndhI", "ndhJ", "ndhK", "psbN", "petA", "petB", "petD", "petG", "petL", "petN", "psaA", "psaB", "psaC", "psaJ", "psaI", "psaM", "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK", "psbL", "psbM", "psbT", "psbZ", "rbcL", "rpl14", "rpl16", "rpl20", "rpl21", "rpl22", "rpl23", "rpl2", "rpl32", "rpl33", "rpl36", "rpoA", "rpoB", "rpoC1", "rpoC2", "rps11", "rps12", "rps14", "rps15", "rps18", "rps19", "rps2", "rps3", "rps4", "rps7", "rps8", "ycf12", "ycf1", "ycf2", "ycf3", "ycf4", "ycf94"]
#locus_list_mt = ["atp1", "atp4", "atp6", "atp8", "atp9", "cob", "cox1", "cox2", "cox3", "matR", "mttB", "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6", "nad7", "nad9", "rpl16", "rpl2", "rpl5", "rpl6", "rps10", "rps11", "rps12", "rps13", "rps14", "rps19", "rps1", "rps2", "rps3", "rps4", "rps7", "sdh3", "sdh4"]
# coding_list = ["rbcL", "rps4"] 
# noncoding_list = ["rps4_trnS", "trnLLF"]
loci_list = coding_list# + noncoding_list

for locus in loci_list:
	locus_file_name = import_align_dir + locus + "_aligned.fasta" #file type
	#print(locus)
	loci_taxa_presence_dict[locus] = ""
	fa_ID = []
	fa_Seq = []
	fa_Num = -1
	for line in open(locus_file_name,"r").readlines():
		line = line.rstrip()
		if line.startswith('>'):
			fa_ID.append(line.replace('>',''))
			fa_Num = fa_Num + 1
			fa_Seq.append("")
		else:
			fa_Seq[fa_Num] = fa_Seq[fa_Num] + line
	fastaID_seq_dic = dict(zip(fa_ID,fa_Seq))
	alignment_len = len(fa_Seq[0])
	for taxon in taxa_list:
		try: 
			taxa_to_fastaID_seq_dic[taxon] = taxa_to_fastaID_seq_dic[taxon] + str(fastaID_seq_dic[taxon])
			loci_taxa_presence_dict[locus] = loci_taxa_presence_dict[locus] + "p\t"
		except:
			taxa_to_fastaID_seq_dic[taxon] = taxa_to_fastaID_seq_dic[taxon] + "-"*alignment_len
			loci_taxa_presence_dict[locus] = loci_taxa_presence_dict[locus] + "ab\t"
	#print locus + ' = ' + str(cumulative_len + 1) + '-' + str(cumulative_len + alignment_len) + ';'
	if locus in coding_list:
		CODON1 = 'DNA, ' + locus + '_pos1 = ' + str(cumulative_len + 1) + '-' + str(cumulative_len + alignment_len) + '\\3'
		CODON2 = 'DNA, ' + locus + '_pos2 = ' + str(cumulative_len + 2) + '-' + str(cumulative_len + alignment_len) + '\\3'
		CODON3 = 'DNA, ' + locus + '_pos3 = ' + str(cumulative_len + 3) + '-' + str(cumulative_len + alignment_len) + '\\3'
		output_parition_file.write(CODON1 + "\n" + CODON2 + "\n" + CODON3 + "\n")
		output_parition12_file.write(CODON1 + "\n" + CODON2 + "\n")
		output_parition3_file.write(CODON3 + "\n")
#	if locus in noncoding_list:
#		NONCODE = 'DNA, ' + locus + ' = ' + str(cumulative_len + 1) + '-' + str(cumulative_len + alignment_len)
#		output_parition_file.write(NONCODE + "\n")
	cumulative_len = cumulative_len + alignment_len
output_parition_file.close()
output_parition12_file.close()
output_parition3_file.close()

for taxon in sorted(taxa_to_fastaID_seq_dic):
	output_alignment.write(">" + taxon + "\n" + taxa_to_fastaID_seq_dic[taxon] + "\n")
output_alignment.close

output_presence_table.write('\t')
for taxon in taxa_list:
	output_presence_table.write(taxon + '\t')
output_presence_table.write('\n')

for locus in sorted(loci_taxa_presence_dict):
	output_presence_table.write(locus + '\t' + loci_taxa_presence_dict[locus] + '\n')


#print taxa_to_fastaID_seq_dic

