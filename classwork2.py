from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import argparse

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Fasta')
    
    parser.add_argument('-i', '--infile', help='input fasta file', type=str, required=True)
    parser.add_argument('-o', '--outfile', help='output fasta file', type=str, required=True)
    
    
    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile


    with open(infile, 'r') as handle:
        alignments_dict = {}
        for record in SeqIO.parse(handle, 'fasta'):
        
            read = str(record.seq)
            results = NCBIWWW.qblast('blastn', 'nt', read, format_type='XML', alignments=1, descriptions=1)
            parser = NCBIXML.parse(results)
            
            for rec in parser:
                if rec.alignments[0].title.split('|')[4] in alignments_dict:
                    alignments_dict[rec.alignments[0].title.split('|')[4]].append(read)
                else:
                    alignments_dict[rec.alignments[0].title.split('|')[4]] = [read]
                

    with open(outfile, 'w') as out:
        for key in alignments_dict.keys():
            if len(alignments_dict[key]) > 1:
                a = 0
                for seq in alignments_dict[key]:
                    a += 1
                    out.write('>' + key + ' ' + str(a) + '\n')
                    out.write(seq + '\n')
            else:
                out.write('>' + key + '\n')
                out.write(alignments_dict[key][0] + '\n') 
    
