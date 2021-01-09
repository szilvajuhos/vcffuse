from vcffuse import ExonCoords
from intervaltree import Interval, IntervalTree


class FusionMaker:
    """
    Joins coordinates of exons for two genes
    """
    # directions from the breakpoint
    DIR_LEFT = True
    DIR_RIGHT = False

    def __init__(self, p5, p3, start, end, chrom):
        """
        :param chrom: chromosome of the call (first col in the VCF)
        """
        self.prime5 = ExonCoords.ExonCoords(p5.chromosome, p5.strand, 0, p5.gene_name, p5.exons)
        self.prime3 = ExonCoords.ExonCoords(p3.chromosome, p3.strand, 0, p3.gene_name, p3.exons)

        # now assign breakpoints to genes:
        # since using the Manta VCF line it is not yet clear
        # which gene contain which endpoint, we have to assign
        # them separately - only if they are meaningful values
        if start is not None and end is not None and start > 0 and end > 0:
            self.assign_breakpoint(start, chrom)
            self.assign_breakpoint(end, chrom)

    def assign_breakpoint(self, bp, chrom):
        for gene in [self.prime5, self.prime3]:
            gene_coords = gene.exons
            gene_ends = IntervalTree()
            gene_ends.add(Interval(gene_coords.begin(), gene_coords.end()))
            if gene_ends.at(bp):
                gene.breakpoint = bp
                print("Breakpoint", bp, "assigned to", gene.gene_name)
                break;
            else:
                print("Breakpoint", bp, "is outside gene", gene.gene_name)
        # there can be a chance that the breakpoint is outside of both genes, assign it to the closest one
        for gene in [self.prime5, self.prime3]:
            if gene.breakpoint == 0 and gene.chromosome == chrom and (
                    abs(gene.exons.begin() - bp) < 20000 or abs(gene.exons.end() - bp) < 20000):
                gene.breakpoint = bp
                print("Breakpoint", bp, "assigned to", gene.gene_name)

    def assign_breakpoint_to_genes(self, bp: tuple):
        chromosome = bp[0]
        for gene in [self.prime5, self.prime3]:
            if chromosome == gene.chromosome:
                gene_extremes = IntervalTree()
                gene_extremes.add(Interval(gene.exons.begin(), gene.exons.end()))
                if gene_extremes.at(bp[1]):
                    gene.breakpoint = bp[1]

    def get_left_part(self, gene: ExonCoords):
        # |-->---!->-->-->--|
        # xxxxxxxx
        # if the breakpoint is in an intron, we have to add a 1-base long interval at the breakpoint
        breakpoint_in_exon = True if len(gene.exons.at(gene.breakpoint)) > 0 else False
        if not breakpoint_in_exon:
            print("**** intron breakpoint at -> ", gene.chromosome + ":" + str(gene.breakpoint - 1))
            gene.exons.add(Interval(gene.breakpoint - 1, gene.breakpoint))
        # when the breakpoint is in an exon, we have to shorten that one
        # by chopping out the unneeded part
        # TODO check it for tandem repeats as looks it was a bug
        # it is still a bug actually
        gene.exons.chop(gene.breakpoint, gene.exons.end())

        return gene.exons

    def get_right_part(self, gene: ExonCoords):
        # |--<---!-<--<--<--|
        #        xxxxxxxxxxxx
        tr_exons = gene.exons
        tr_exons = tr_exons.overlap(gene.breakpoint, gene.exons.end())
        tr_exons.add(Interval(gene.breakpoint, gene.breakpoint + 1))
        return tr_exons

    def print_as_bed(self, chromosome, exs):
        for e in sorted(exs):
            print(chromosome + "\t" + str(e.begin) + "\t" + str(e.end))

    def getFusedPart(self, gene: ExonCoords, direction) -> ExonCoords:
        """
        Breaks the gene coordinates at the breakpoint
        and returs with the left or right truncated part
        """
        tr_exons = None
        if direction == 5:
            if gene.strand > 0:
                tr_exons = self.get_left_part(gene)
            else:
                tr_exons = self.get_right_part(gene)
        else:  # assuming prime 3 - the other way around
            if gene.strand < 0:
                tr_exons = self.get_right_part(gene)
            else:
                tr_exons = self.get_left_part(gene)
        return ExonCoords(gene.chromosome, gene.strand, gene.breakpoint, gene.gene_name, tr_exons)

    def fuse_deletion(self):
        prime5part = None
        prime3part = None
        if self.prime5.strand < 0:
            # |<-<-<-<-|     |<-<-<-<-|
            #     |-------------|
            # we need right part for 5' and left part for 3'
            prime5part = ExonCoords(self.prime5.chromosome,
                                    self.prime5.strand,
                                    self.prime5.breakpoint,
                                    self.prime5.gene_name,
                                    self.get_right_part(self.prime5))
            prime3part = ExonCoords(self.prime3.chromosome,
                                    self.prime3.strand,
                                    self.prime3.breakpoint,
                                    self.prime3.gene_name,
                                    self.get_left_part(self.prime3))
        else:
            # |->->->->|     |->->->->|
            #     |-------------|
            # we need left part for 5' and right part for 3'
            prime5part = ExonCoords(self.prime5.chromosome,
                                    self.prime5.strand,
                                    self.prime5.breakpoint,
                                    self.prime5.gene_name,
                                    self.get_left_part(self.prime5))
            prime3part = ExonCoords(self.prime3.chromosome,
                                    self.prime3.strand,
                                    self.prime3.breakpoint,
                                    self.prime3.gene_name,
                                    self.get_right_part(self.prime3))
        prime5part.print_as_bed()
        prime3part.print_as_bed()
        # now move the 3' part to the 5' part
        p5borders = (prime5part.exons.begin(), prime5part.exons.end())
        p3borders = (prime3part.exons.begin(), prime3part.exons.end())
        shift = 0
        if prime5part.strand > 0:
            shift = prime5part.breakpoint - prime3part.exons.begin() + 1
        else:
            shift = prime5part.exons.begin() - prime3part.exons.end()
        # we have to shift 3' only
        shifted3p = IntervalTree()
        for iv in prime3part.exons:
            shifted3p.add(Interval(iv.begin + shift, iv.end + shift))
        shifted5p = prime5part.exons
        # and now shift down the stuff to 0 for SVG
        left_shift = (shifted5p | shifted3p).begin()
        # TODO: DRY it out
        based05p = IntervalTree()
        for iv in shifted5p:
            based05p.add(Interval(iv.begin - left_shift, iv.end - left_shift))
        based03p = IntervalTree()
        for iv in shifted3p:
            based03p.add(Interval(iv.begin - left_shift, iv.end - left_shift))
        return based05p, based03p

    def fuse_tandem_genes(self):
        # dealing with tandem repeats:
        # first we have to get parts by strand
        # - strand means we want to have the
        #       right part from the breakpoint for 5'
        #       left part from the breakpoint for 3'
        #       join them by starting with the 5' part,
        #       add the 3' part to its left
        # + strand means we want to have the
        #       left part from the breakpoint for 5'
        #       right part from the breakpoint for 3'
        #       join them by starting with the 5' part
        #       add the 3' part to its right
        prime5part = None
        prime3part = None
        if (self.prime5.strand < 0):
            prime5part = ExonCoords(self.prime5.chromosome,
                                    self.prime5.strand,
                                    self.prime5.breakpoint,
                                    self.prime5.gene_name,
                                    self.get_right_part(self.prime5))
            prime3part = ExonCoords(self.prime3.chromosome,
                                    self.prime3.strand,
                                    self.prime3.breakpoint,
                                    self.prime3.gene_name,
                                    self.get_left_part(self.prime3))
        else:
            prime5part = ExonCoords(self.prime5.chromosome,
                                    self.prime5.strand,
                                    self.prime5.breakpoint,
                                    self.prime5.gene_name,
                                    self.get_left_part(self.prime5))
            prime3part = ExonCoords(self.prime3.chromosome,
                                    self.prime3.strand,
                                    self.prime3.breakpoint,
                                    self.prime3.gene_name,
                                    self.get_right_part(self.prime3))
        prime5part.print_as_bed()
        prime3part.print_as_bed()
        # now move the 3' part to the 5' part
        p5borders = (prime5part.exons.begin(), prime5part.exons.end())
        p3borders = (prime3part.exons.begin(), prime3part.exons.end())
        # |------5------|
        #                   |------3------|
        # |------5------|
        #           |------3------|
        #
        #                   |------5------|
        #           |------3------|
        #                   |------5------|
        # |------3------|
        # shift = (5'start - 3'start)
        # 3'start = 3'start + shift
        shift = 0
        if prime5part.strand > 0:
            shift = prime5part.breakpoint - prime3part.exons.begin() + 1
        else:
            shift = prime5part.exons.begin() - prime3part.exons.end()
        # we have to shift 3' only
        shifted3p = IntervalTree()
        for iv in prime3part.exons:
            shifted3p.add(Interval(iv.begin + shift, iv.end + shift))
        shifted5p = prime5part.exons
        # and now shift down the stuff to 0 for SVG
        left_shift = (shifted5p | shifted3p).begin()
        # TODO: DRY it out
        based05p = IntervalTree()
        for iv in shifted5p:
            based05p.add(Interval(iv.begin - left_shift, iv.end - left_shift))
        based03p = IntervalTree()
        for iv in shifted3p:
            based03p.add(Interval(iv.begin - left_shift, iv.end - left_shift))
        return based05p, based03p

    def fuse_inversion(self):
        based05p = None
        based03p = None
        if self.prime5.strand > 0:  # forward 5'
            print("Forward 5' inversion")
            prime5part = ExonCoords(self.prime5.chromosome, self.prime5.strand,
                                    self.prime5.breakpoint, self.prime5.gene_name,
                                    self.get_left_part(self.prime5))
            prime3part = ExonCoords(self.prime3.chromosome, self.prime3.strand,
                                    self.prime3.breakpoint, self.prime3.gene_name,
                                    self.get_left_part(self.prime3))
            # we have to turn around the 3' part
            # -->  <--
            prime3part = self.turn_backwards(prime3part)
            # -->  -->
        else:
            print("Reverse 5' inversion")
            prime5part = ExonCoords(self.prime5.chromosome, self.prime5.strand,
                                    self.prime5.breakpoint, self.prime5.gene_name,
                                    self.get_right_part(self.prime5))
            prime3part = ExonCoords(self.prime3.chromosome, self.prime3.strand,
                                    self.prime3.breakpoint, self.prime3.gene_name,
                                    self.get_right_part(self.prime3))
            # now the inversion is like
            # <--  -->
            prime5part = self.turn_backwards(prime5part)
            # -->  -->
        based05p = self.shift_left_to(0, prime5part)
        based03p = self.shift_left_to(based05p.end(), prime3part)
        print(based05p)
        print(based03p)
        return based05p, based03p

    def fuse_translocations(self, p5dir, p3dir):
        prime5part = None
        prime3part = None
        # TODO: DRY it out
        if p5dir == self.DIR_LEFT:
            prime5part = ExonCoords(self.prime5.chromosome, self.prime5.strand,
                                    self.prime5.breakpoint, self.prime5.gene_name,
                                    self.get_left_part(self.prime5))
        else:
            prime3part = ExonCoords(self.prime3.chromosome, self.prime3.strand,
                                    self.prime3.breakpoint, self.prime3.gene_name,
                                    self.get_right_part(self.prime3))
        if p3dir == self.DIR_LEFT:
            prime3part = ExonCoords(self.prime3.chromosome, self.prime3.strand,
                                    self.prime3.breakpoint, self.prime3.gene_name,
                                    self.get_left_part(self.prime3))
        else:
            prime3part = ExonCoords(self.prime3.chromosome, self.prime3.strand,
                                    self.prime3.breakpoint, self.prime3.gene_name,
                                    self.get_right_part(self.prime3))

        if p5dir == self.DIR_LEFT and p3dir == self.DIR_LEFT:
            # forward antiparallel
            # we have to turn the reverse strand 3' gene backwards
            prime3part = self.turn_backwards(prime3part)
            # and have to stick it to the 5' part
            # we are ignoring chromosomes this time
            # TODO: DRY it out
            based05p = self.shift_left_to(0, prime5part)
            based03p = self.shift_left_to(based05p.end(), prime3part)
            print(based05p)
            print(based03p)
            return (based05p, based03p)

    def shift_left_to(self, new_base: int, exs: ExonCoords):
        rebased = IntervalTree()
        shift = exs.exons.begin() - new_base
        for item in exs.exons:
            rebased.add(Interval(item.begin - shift, item.end - shift))
        return rebased

    def turn_backwards(self, exs: ExonCoords):
        """
        We have coords like:
        |###|----|####|--|#|
        and want to have something like:
        |#|--|####|----|###|
        :param exs:
        ExonCoords that we want to turn backwards
        :return:
        new IntervalTree() with backwards coordinates
        """
        gene_start = exs.exons.begin()
        gene_end = exs.exons.end()
        new_exons = IntervalTree()
        for item in sorted(exs.exons):
            new_end = gene_end - (item.begin - gene_start)
            new_start = gene_end + gene_start - item.end
            # print(exs.chromosome + "\t" + str(new_start) + "\t" + str(new_end) + "\t"+ exs.gene_name)
            new_exons.add(Interval(new_start, new_end))
        return ExonCoords(exs.chromosome, exs.strand, exs.breakpoint, exs.gene_name, new_exons)

    def print_properties(self):
        print("5' gene:")
        print("strand     :", self.prime5.strand)
        print("breakpoint :", self.prime5.chromosome + ":" + str(self.prime5.breakpoint))
        exons = self.prime5.exons
        start = exons.begin()
        end = exons.end()
        print("gene coords:", self.prime5.chromosome + ":" + str(start) + "-" + str(end))
        print("3' gene:")
        print("strand     :", self.prime3.strand)
        print("breakpoint :", self.prime3.chromosome + ":" + str(self.prime3.breakpoint))
        exons = self.prime3.exons
        start = exons.begin()
        end = exons.end()
        print("gene coords:", self.prime3.chromosome + ":" + str(start) + "-" + str(end))
        self.prime5.print_as_bed()
        self.prime3.print_as_bed()
