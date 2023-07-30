from unittest import TestCase
import pytest
import os, shutil, tempfile

from mapper import Read, Reference, Mapping, read_fasta, SAMWriter, ReadPolisher


class MapperTestCase(TestCase):
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    @pytest.mark.read
    def test_read_constructor(self):
        read = Read([">Read_0", "AGTCGTAG", "TTCAGCCT", "CGTTAGCT", "AGGCAATG"])
        self.assertEqual("AGTCGTAGTTCAGCCTCGTTAGCTAGGCAATG", read.bases)
        self.assertEqual("Read_0", read.name)

    @pytest.mark.read
    def test_read_str(self):
        read = Read([">Read_0", "AGTCGTAG", "TTCAGCCT", "CGTTAGCT", "AGGCAATG"])
        self.assertEqual("Read_0: AGTCGTAGTTCAGCCTCGTT...", str(read))

    @pytest.mark.read
    def test_read_seed(self):
        read = Read([">Read_0", "AGTCGTAG", "TTCAGCCT", "CGTTAGCT", "AGGCAATG"])
        self.assertEqual("AGTCG", read.get_seed(5))

    @pytest.mark.reference
    def test_reference_constructor(self):
        read = Reference([">Reference", "TTTACTGTGTCCATGGTGTATCCTGTTCCT", "GTTCCATGGCTGTATGGAGGATCTCCAGTATAAGAGAATG"])
        self.assertEqual("TTTACTGTGTCCATGGTGTATCCTGTTCCTGTTCCATGGCTGTATGGAGGATCTCCAGTATAAGAGAATG", read.bases)
        self.assertEqual("Reference", read.name)

    @pytest.mark.reference
    def test_reference_str(self):
        read = Reference([">Reference", "TTTACTGTGTCCATGGTGTATCCTGTTCCT", "GTTCCATGGCTGTATGGAGGATCTCCAGTATAAGAGAATG"])
        self.assertEqual("Reference: TTTACTGTGTCCATGGTGTA...", str(read))

    @pytest.mark.reference
    def test_reference_get_kmer_positions(self):
        ref = Reference([">ref", "AGTCCTGATTAGCGGTTAGCGAAT"])
        self.assertEqual([9, 16], ref.get_kmer_positions("TAG"))
        self.assertEqual([0], ref.get_kmer_positions("AGTC"))

    @pytest.mark.reference
    def test_reference_count_mismatches(self):
        ref = Reference([">ref", "AGTCCTGATTAGCGGTTAGCGAAT"])
        read1 = Read([">read_1", "CCTGAT"])
        self.assertEqual(4, ref.count_mismatches(read1, 0))
        self.assertEqual(0, ref.count_mismatches(read1, 3))

    @pytest.mark.toplevel
    def test_read_fasta_reference(self):
        references = read_fasta("data/fluA.fasta", Reference.__name__)
        self.assertEqual(1, len(references))
        reference = references[0]
        self.assertEqual(Reference.__name__, reference.__class__.__name__)
        self.assertEqual(len(reference.bases), 2445)

    @pytest.mark.toplevel
    def test_read_fasta_reads(self):
        reads = read_fasta("data/fluA_reads.fasta", Read.__name__)
        self.assertEqual(100, len(reads))
        read = reads[12]
        self.assertEqual(Read.__name__, read.__class__.__name__)
        self.assertEqual(read.bases, "GGCAAAAATAATGAATTTAACTTGTCCTTCATGAAAAAATGCCTGTTTTT")

    @pytest.mark.mapping
    def test_mapping(self):
        ref = Reference([">ref", "AGTCCTGATTAGCGGTTAGCGAAT"])
        read1 = Read([">read_1", "CCTGAT"])
        read2 = Read([">read_2", "TAGCGGT"])
        mapping = Mapping(ref)
        self.assertEqual([], mapping.get_reads_at_position(0))
        mapping.add_read(read1, 5)
        self.assertEqual([read1], mapping.get_reads_at_position(5))
        self.assertEqual([], mapping.get_reads_at_position(0))
        mapping.add_read(read2, 5)
        self.assertEqual([read1, read2], mapping.get_reads_at_position(5))

    @pytest.mark.samwriter
    def test_writer(self):
        ref = Reference([">ref additional reference information", "AGTCCTGATTAGCGGTTAGCGAAT"])
        read1 = Read([">read_1", "CCTGAT"])
        read2 = Read([">read_2", "TAGCGGT"])
        mapping = Mapping(ref)
        mapping.add_read(read1, 5)
        mapping.add_read(read2, 9)
        w = SAMWriter(mapping)
        outname = os.path.join(self.test_dir, "out.sam")
        w.write_mapping(outname)
        with open(outname, "r") as f:
            lines = []
            for line in f.readlines():
                lines += [line.strip()]
        assert lines[0].strip()=="@SQ\tSN:ref\tLN:24"
        assert "read_1\t0\tref\t6\t255\t6M\t*\t0\t0\tCCTGAT\t*" in lines
        assert "read_2\t0\tref\t10\t255\t7M\t*\t0\t0\tTAGCGGT\t*" in lines
        assert len([x for x in lines if len(x) != 0]) == 3

    @pytest.mark.polisher
    def test_replacements(self):
        p = ReadPolisher(3)
        p.add_read("AGTCG")
        p.add_read("AGTCG")
        p.add_read("ACTCG")
        rep = p.get_replacements(3)
        assert rep == {}
        rep = p.get_replacements(2)
        assert rep == {'CTC': 'GTC', 'ACT': 'AGT'}

    @pytest.mark.polisher
    def test_replace(self):
        r = Read([">r1", "ACTCG"])
        r.replace_kmers({'CTC': 'GTC', 'ACT': 'AGT'})
        assert r.get_seed(3) == "AGT"