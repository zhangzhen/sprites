BASEPAIRING = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'g':'c', 'c':'g'}

class Seq(object):
    def __init__(self, name, text):
        self.name = name
        self.text = text

class Event(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
    def __str__(self):
        return "%s\t%d\t%d" % (self.chrom, self.start, self.end)

class Deletion(Event):
    def __str__(self):
        return Event.__str__(self) + "\tDEL\t."

class Inversion(Event):
    def __str__(self):
        return Event.__str__(self) + "\tINV\t."

class Insertion(Event):
    def __init__(self, chrom, start, end, text):
        Event.__init__(self, chrom, start, end)
        self.text = text
    def __str__(self):
        return Event.__str__(self) + "\tINS\t%s" % self.text
        
def reversecomplement(seq):
    comp = [BASEPAIRING[ch] for ch in seq]
    comp.reverse()
    return "".join(comp)

def buildchrom(ref, events):
    pos = 0
    chunks = []
    text = ref.text
    for e in events:
        chunks.append(text[pos:e.start-1])
        if isinstance(e, Insertion):
            chunks.append(e.text)
        if isinstance(e, Inversion):
            chunks.append(reversecomplement(text[e.start-1:e.end-1]))
        pos = e.end - 1
    chunks.append(text[pos:])
    return Seq(ref.name, "".join(chunks))
        
def tofile(out, events):
    del_cnt = ins_cnt = inv_cnt = 0
    out.write("chr\tstart\tend\tclass\tseq\n")
    for e in events:
        if isinstance(e, Deletion):
            del_cnt += 1
        if isinstance(e, Insertion):
            ins_cnt += 1
        if isinstance(e, Inversion):
            inv_cnt += 1
        out.write("%s\n" % str(e))
    out.close()
    print "\nSummary:\n# %d deletions, %d insertions and %d inversions in the generated sequence\n" \
              % (del_cnt, ins_cnt, inv_cnt)

def loadevents(reader, chrom):
    events = []
    for e in reader:
        if e is None or chrom != e.chrom:
            continue
        events.append(e)
    return events

def isgffsorted(reader):
    prev = reader.next()
    for e in reader:
        if e.chrom < prev.chrom:
            return False
        if e.chrom == prev.chrom and e.start < prev.start:
            return False
        prev = e
    return True

class GffReader(object):
    def __init__(self, file):
        self.file = file

    def __iter__(self):
        return self

    def next(self):
        line = self.file.readline()
        if '' == line:
            self.file.close()
            raise StopIteration
        r = line.rstrip('\n').split('\t')
        if r[2] == "homozygous_indel" and r[10] == "Homozygous_Deletion":
            return Deletion(r[0], int(r[3]), int(r[4]))
        if r[2] == "homozygous_indel" and r[10] == "Homozygous_Insertion":
            s = reversecomplement(r[9]) if r[6] == '-' else r[9]
            return Insertion(r[0], int(r[3]), int(r[4]), s)
        if r[2] == "assembly_comparison_inversion":
            return Inversion(r[0], int(r[3]), int(r[4]))
        return None
