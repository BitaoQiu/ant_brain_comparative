
#1. Transcript2gene, could make a dict and add prefix.
#2. Replace the ortholog table with ortholog gene table.

def get_info(info): #Using gff file to format the transcript and gene relation.
    info = info.split(';')
    transcript_ID = [i for i in info if i.startswith('ID=')][0].replace('ID=','').rstrip('\n').upper()
    gene_ID = [i for i in info if i.startswith('Parent=')][0].replace('Parent=','').rstrip('\n')
    return transcript_ID,gene_ID

def gff2dict(prefix, sp_dict):
    sp_dict[prefix] = {}
    gff_input = open(prefix+'.gff','r').readlines()
    for line in gff_input[1:]:
        line = line.split('\t')
        if line[2] == 'mRNA':
            transcript_ID,gene_ID = get_info(line[8])
            gene_ID = prefix + "_" + gene_ID
            sp_dict[prefix][transcript_ID] = gene_ID
    return 

t2g_dict = {}
gff2dict('Aech',t2g_dict)
gff2dict('Mpha',t2g_dict)
gff2dict('Lhum',t2g_dict)
gff2dict('Sinv',t2g_dict)
gff2dict('Dqua',t2g_dict)
gff2dict('Cbir',t2g_dict)
gff2dict('Lnig',t2g_dict)

ortholog_table = open(r'myproject.poff','r').readlines()
ortholog_gene_table = open(r'gene_table.poff','w')
def output_gene(line, sp_dict):
    Aech = sp_dict['Aech'][line[3]]
    Cbir = sp_dict['Cbir'][line[4]]
    Dqua = sp_dict['Dqua'][line[5]]
    Lhum = sp_dict['Lhum'][line[6]]
    Lnig = sp_dict['Lnig'][line[7]]
    Mpha = sp_dict['Mpha'][line[8]]
    Sinv = sp_dict['Sinv'][line[9]]
    return Aech, Sinv, Mpha, Lnig, Lhum, Cbir, Dqua
    
for line in ortholog_table[1:]:
    line = line.rstrip('\n').split('\t')
    if line[0] == '7' and line[1] == '7': # Only retained genes wit one-to-one orthologos across seven ant species.
        print('\t'.join(output_gene(line, t2g_dict)),file = ortholog_gene_table)
ortholog_gene_table.close()
    

