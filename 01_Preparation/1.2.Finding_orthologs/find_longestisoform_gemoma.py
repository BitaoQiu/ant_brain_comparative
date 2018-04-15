'''Find the longest isoform protein for each gene, and change the ID to uppercase (because it's upper case in fasta file!!!!'''
import sys
from Bio import SeqIO
gff_file = open(sys.argv[1], 'r').readlines()
pep_file = sys.argv[2]
out_gff = open(sys.argv[3] + '.gff','w')
out_fasta = open(sys.argv[3] + '.fasta','w')
pep_seq_dict = SeqIO.to_dict(SeqIO.parse(pep_file, "fasta"))

def info_grep(info):
    info_all = info.split(';')
    mRNA_ID = [i for i in info_all if i.startswith('ID=')][0][3:].upper()
    parent_ID = [i for i in info_all if i.startswith('Parent=')][0][7:].rstrip('\n')
    AA_length = int([i for i in info_all if i.startswith('AA=')][0][3:])
#    AA_coverage = float([i for i in info_all if i.startswith('avgCov=')][0][7:])
    return mRNA_ID, parent_ID, AA_length#, AA_coverage

gene_dict = dict()

for line in gff_file[1:]:
    line = line.split('\t')
    info = line[8]
    if line[2] == 'mRNA':
        mRNA_ID, parent_ID, AA_length = info_grep(info)#, AA_coverage = info_grep(info)
        if parent_ID not in gene_dict:
            gene_dict[parent_ID] = {}
        if mRNA_ID in pep_seq_dict:
            gene_dict[parent_ID][mRNA_ID] = {}
            gene_dict[parent_ID][mRNA_ID]['length'] = AA_length
#            gene_dict[parent_ID][mRNA_ID]['coverage'] = AA_coverage
            gene_dict[parent_ID][mRNA_ID]['other_info'] = line

for gene in gene_dict:
    try:
        longest_var = sorted(gene_dict[gene], key = lambda x: gene_dict[gene][x]['length'])[-1]
        out_info = gene_dict[gene][longest_var]['other_info']
        out_info[2] = 'CDS'
        out_name = 'ID='+longest_var+';'
        out_info[8] = out_name
        print('\t'.join(out_info), file = out_gff)
        print('>'+longest_var, file = out_fasta)
        print(str(pep_seq_dict[longest_var].seq), file = out_fasta)
    except:
        print(gene)

out_gff.close()
out_fasta.close()


