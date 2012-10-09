#!/usr/bin/python
import sys
import ParseFastQ

def filterReads(inFilename, outFilename, numN=5):
    parser = ParseFastQ.ParseFastQ(inFilename)  # optional arg: headerSymbols allows changing the header symbols
    output = open(outFilename, 'w')
    total = 0
    rest = 0
    for record in parser:
        total += 1
        if 'N'*numN in record[1]:
            continue
        output.write("%s\n%s\n%s\n%s\n" % (record[0], record[1], record[2], record[3]))
        rest += 1
    print "total #reads: %d" % total
    print "remaining #reads: %d" % rest

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Usage: %s [input_file.fa] [output_file.fa]" % sys.argv[0]
        sys.exit(1)
    filterReads(sys.argv[1], sys.argv[2])

        
    

