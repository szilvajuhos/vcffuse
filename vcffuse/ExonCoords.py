from intervaltree import Interval, IntervalTree


class ExonCoords:
    def __init__(self, chromosome, strand, breakpoint, gene_name, exons: IntervalTree):
        self.chromosome = chromosome
        self.strand = strand
        self.breakpoint = breakpoint
        self.gene_name = gene_name
        self.exons = IntervalTree(exons)

    @classmethod
    def fromTuple(cls, a_tuple):
        return cls(a_tuple[0], a_tuple[1], a_tuple[2], a_tuple[3], a_tuple[4])

    @classmethod
    def copy_without_exons(cls, exc):
        return cls(exc.chromosome, exc.strand, exc.breakpoint, exc.gene_name, IntervalTree())

    @classmethod
    def empty(cls):
        return cls("", 0, -1, "", IntervalTree())

    def print_properties(self):
        print("#########################################")
        print("coordinates :", self.chromosome + ":" + str(self.exons.begin()) + "-" + str(self.exons.end()))
        print("gene        :", self.gene_name)
        print("strand      :", self._strand)
        print("breakpoint  :", self._breakpoint)
        print("exons       :", self._exons)
        print("#########################################")

    def print_as_bed(self):
        chromosome = self.chromosome
        for ex in sorted(self.exons):
            print(chromosome + "\t" + str(ex.begin) + "\t" + str(ex.end))

    @property
    def gene_name(self):
        return self._gene_name

    @gene_name.setter
    def gene_name(self, value):
        self._gene_name = value

    @property
    def chromosome(self):
        return self._chromosome

    @chromosome.setter
    def chromosome(self, value):
        self._chromosome = value

    @property
    def strand(self):
        return self._strand

    @strand.setter
    def strand(self, value):
        self._strand = value

    @property
    def breakpoint(self):  # int
        return self._breakpoint

    @breakpoint.setter
    def breakpoint(self, value):
        self._breakpoint = value

    @property
    def exons(self):  # IntervalTree()
        return self._exons

    @exons.setter
    def exons(self, exons):
        self._exons = exons

    def begin(self):
        return self.exons.begin()
