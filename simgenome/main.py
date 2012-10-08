from bx.seq.fasta import FastaFile, FastaWriter
import sv
import argparse

def main():
    parser = argparse.ArgumentParser(description='Sequence simulator')
    parser.add_argument('-i','--input', help='Sequence fasta file', required=True)
    parser.add_argument('-sv','--svfile', help='Structural variations that is to be inserted', required=True)
    parser.add_argument('-o','--output', help='Sequence with structural variation', required=True)
    parser.add_argument('-r','--report', help='List of structural variations', required=True)
    parser.add_argument('-c','--chrom', help='Chromosome name', required=True)
    args = vars(parser.parse_args())    
    events = sv.loadevents(sv.GffReader(file(args['svfile'], "r")), args['chrom'])
    sv.tofile(file(args['report'], "w"), events)
    ref = FastaFile(file(args['input'], "r"))
    seq = sv.buildchrom(ref, events)
    writer = FastaWriter(file(args['output'], "w"))
    writer.write(seq)

if __name__ == '__main__':
    main()
