from bx.seq.fasta import FastaFile, FastaWriter
import sv
import argparse

def main():
    parser = argparse.ArgumentParser(description='Description of your program')
    parser.add_argument('-i','--input', help='Description for foo argument', required=True)
    parser.add_argument('-sv','--svfile', help='Description for bar argument', required=True)
    parser.add_argument('-o','--output', help='Description for foo argument', required=True)
    parser.add_argument('-r','--report', help='Description for foo argument', required=True)
    parser.add_argument('-c','--chrom', help='Description for foo argument', required=True)
    args = vars(parser.parse_args())    
    events = sv.loadevents(sv.GffReader(file(args['svfile'], "r")), args['chrom'])
    sv.tofile(file(args['report'], "w"), events)
    ref = FastaFile(file(args['input'], "r"))
    seq = sv.buildchrom(ref, events)
    writer = FastaWriter(file(args['output'], "w"))
    writer.write(seq)

if __name__ == '__main__':
    main()
