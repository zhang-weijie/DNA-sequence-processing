import os


class Sequence:
    def __init__(self, lines):
        self.name = lines[0].strip()[1:]
        self.bases = "".join([x.strip() for x in lines[1:]]).upper()

    def __str__(self):
        return self.name + ": " + self.bases[:20] + "..."

    def __repr__(self):
        return self.__str__()


class Read(Sequence):
    def get_seed(self, seedlength):
        return self.bases[:seedlength]

    def replace_kmers(self, replacements):
        for key in replacements:
            self.bases = self.bases.replace(key, replacements[key])

class Reference(Sequence):
    def __init__(self, lines):
        self.kmers = None
        super().__init__(lines)

    def calculate_kmers(self, kmersize):
        self.kmers = {}
        for pos in range(0, len(self.bases) - kmersize + 1):
            kmer = self.bases[pos:(pos + kmersize)]
            if kmer not in self.kmers:
                self.kmers[kmer] = []
            self.kmers[kmer] += [pos]

    def get_kmer_positions(self, kmer):
        if self.kmers is None or len(next(iter(self.kmers))) != len(kmer):
            self.calculate_kmers(len(kmer))
        if kmer not in self.kmers:
            return []
        return self.kmers[kmer]

    def count_mismatches(self, read, position):
        mismatches = 0
        for pos in range(position, position + len(read.bases)):
            if pos >= len(self.bases):
                break
            if read.bases[pos - position] != self.bases[pos]:
                mismatches += 1
        # Count every base of the read that goes out past the end of the reference as a mismatch
        mismatches += position + len(read.bases) - pos - 1
        return mismatches


class Mapping:
    def __init__(self, reference):
        self.reference = reference
        self.reads = {}

    def add_read(self, read, position):
        if position not in self.reads:
            self.reads[position] = []
        self.reads[position] += [read]

    def get_reads_at_position(self, position):
        if position not in self.reads:
            return []
        return self.reads[position]

    def __str__(self):
        res = ["Mapping to " + self.reference.name]
        for pos in self.reads:
            res += ["  " + str(len(self.reads[pos])) + " reads mapping at " + str(pos)]
        return "\n".join(res)


class SAMWriter:
    def __init__(self, mapping):
        self.mapping = mapping

    def write_mapping(self, filename):
        with open(filename, "w") as f:
            ref = self.mapping.reference
            ref_name = ref.name.split(" ")[0]
            ref_length = len(ref.bases)
            header = "".join(["@SQ\tSN:", ref_name, "\tLN:", str(ref_length), "\n"])
            f.write(header)

            reads = self.mapping.reads
            for pos in reads:
                for read in reads[pos]:
                    read_name = read.name
                    read_bases = read.bases
                    cigar = str(len(read_bases)) + "M"
                    sam = "".join(
                        [read_name, "\t0\t", ref_name, "\t", str(pos + 1), "\t255\t", cigar, "\t*\t0\t0\t", read_bases,
                         "\t*\n"])
                    f.write(sam)


class ReadPolisher:
    def __init__(self, kmerlen):
        self.kmerlen = kmerlen
        self.spectrum = {}

    def add_read(self, readseq):
        if len(readseq) < self.kmerlen:
            return
        for i in range(len(readseq) - self.kmerlen + 1):
            subseq = readseq[i: i + self.kmerlen]
            if subseq in self.spectrum:
                self.spectrum[subseq] += 1
            else:
                self.spectrum[subseq] = 1

    def get_replacements(self, minfreq):
        replacements = {}
        for subseq in self.spectrum:
            for i in range(self.kmerlen):
                maxfreq_candidate = (-1, "")
                for base in ["A", "G", "T", "C"]:
                    if base != subseq[i]:
                        candidate = subseq[:i] + base + subseq[i+1:]
                        if candidate in self.spectrum and self.spectrum[candidate] >= minfreq:
                            if self.spectrum[candidate] > maxfreq_candidate[0]:
                                maxfreq_candidate = (self.spectrum[candidate], candidate)
                if maxfreq_candidate[0] != -1:
                    replacements[subseq] = maxfreq_candidate[1]
        return replacements


def read_fasta(fastafile, klassname):
    klass = globals()[klassname]
    f = open(fastafile, "r")
    readlines = []
    reads = []
    for line in f:
        if line[0] == '>' and len(readlines) != 0:
            reads += [klass(readlines)]
            readlines = []
        readlines += [line]
    reads += [klass(readlines)]
    f.close()
    return reads


def map_reads(reads, reference, kmersize, max_mismatches):
    mapping = Mapping(reference)
    reference.calculate_kmers(kmersize)
    for read in reads:
        seed = read.get_seed(kmersize)
        seed_positions = reference.get_kmer_positions(seed)
        for position in seed_positions:
            mismatches = reference.count_mismatches(read, position)
            if mismatches < max_mismatches:
                mapping.add_read(read, position)
    return mapping


def main():
    # reads = read_fasta("data/fluA_reads.fasta", Read.__name__)
    # reference = read_fasta("data/fluA.fasta", Reference.__name__)[0]
    # mapping = map_reads(reads, reference, 8, 5)
    # print("Mapping reads: " + str(len(mapping.reads)))

    reads = read_fasta("data/patient2.fasta", Read.__name__)
    reference = read_fasta("data/rpoB.fasta", Reference.__name__)[0]

    # mapping = map_reads(reads, reference, 30, 15)  # p1 null
    # mapping = map_reads(reads, reference, 150, 15)  # p1 alternative null

    # mapping = map_reads(reads, reference, 15, 5)  # p2 keine Mutation bei 1862
    # mapping = map_reads(reads, reference, 30, 15)  # p2 16 * C1862A und 1 Read ohne Mutation bei 1862
    # mapping = map_reads(reads, reference, 150, 15)  # p2 alternative 4 * C1862A und 1 Read ohne Mutation bei 1862

    # mapping = map_reads(reads, reference, 20, 11)  # p3 20+ * C1402A
    # mapping = map_reads(reads, reference, 150, 15)  # p3 alternative 20+ * C1402A

    # mapping = map_reads(reads, reference, 20, 11)  # p4 1 * T2858G
    # mapping = map_reads(reads, reference, 150, 15)  # p4 alternative 1 * T2858G
    # writer = SAMWriter(mapping)
    # writer.write_mapping("data/patient2_mapping.sam")


    # polisher = ReadPolisher(30) #p1
    # polisher = ReadPolisher(150) #p1 alternative

    polisher = ReadPolisher(15) #p2
    # polisher = ReadPolisher(150) #p2 alternative

    # polisher = ReadPolisher(20) #p3
    # polisher = ReadPolisher(150) #p3 alternative

    # polisher = ReadPolisher(30) #p4
    # polisher = ReadPolisher(150) #p4 alternative

    for read in reads:
        polisher.add_read(read.bases)
    replacements = polisher.get_replacements(3)

    for read in reads:
        read.replace_kmers(replacements)

    # mapping = map_reads(reads, reference, 30, 2)  #p1 null
    # mapping = map_reads(reads, reference, 150, 8)  #p1 alternative null

    mapping = map_reads(reads, reference, 15, 2) #p2 10 * C1862A
    # mapping = map_reads(reads, reference, 150, 8) #p2 alternative 7 * C1862A

    # mapping = map_reads(reads, reference, 20, 3) #p3 20+ * C1402A
    # mapping = map_reads(reads, reference, 150, 8) #p3 alternative 20+ * C1402A

    # mapping = map_reads(reads, reference, 30, 4) #p4 1 * T2858G
    # mapping = map_reads(reads, reference, 150, 15) #p4 alternative 1 * T2858G

    writer = SAMWriter(mapping)
    writer.write_mapping("data/patient2_mapping_corrected.sam")


if __name__ == "__main__":
    main()
