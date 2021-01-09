import svgwrite
from vcffuse import ExonCoords
from svgwrite import cm, mm


class CreateSVGPicture:
    def __init__(self):
        self._width = '100%'
        self._height = '100%'
        self._verbose = False

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        self._verbose = value

    def make_SVG_file(self, fex, svg_file_name, pic_count, prime5name, prime3name):
        outfile = str(pic_count) + "_" + prime5name + "_" + prime3name + "-" + svg_file_name
        dwg = svgwrite.Drawing(filename=outfile, size=(self._width, self._height), debug=True)
        # background
        dwg.add(dwg.rect(insert=(0, 0), size=(self._width, self._height), fill='white', stroke='white'))
        # exons
        shapes = dwg.add(dwg.g(id='shapes', fill='red'))
        # 5' part of exons
        shapes = self.shape_intervals(dwg, shapes, fex[0], 'blue')
        # 3' part of exons
        shapes = self.shape_intervals(dwg, shapes, fex[1], 'red')
        # introns
        shapes.add(dwg.rect(insert=(fex[0].begin() / 100 * mm, 8 * mm),
                            size=(int(abs(fex[0].begin() - fex[0].end())) / 100 * mm, 2 * mm),
                            fill='blue'))
        shapes.add(dwg.rect(insert=(fex[1].begin() / 100 * mm, 8 * mm),
                            size=(int(abs(fex[1].end() - fex[1].begin())) / 100 * mm, 2 * mm),
                            fill='red'))
        dwg.save()
        if self._verbose:
            print("Fusion picture is at", outfile)
        return pic_count + 1

    def shape_intervals(self, dwg, shapes, itv: ExonCoords, color):
        height = 2 * cm
        for iv in itv:
            if abs(iv.end - iv.begin) > 1:
                s = iv.begin / 100
                width = int(abs(iv.end - iv.begin + 100) / 100)
                # print("SVG:", s*mm, 0, width*mm, height)
                shapes.add(dwg.rect(insert=(s * mm, 0), size=(width * mm, height),
                                    fill=color, stroke=color, stroke_width=1))
        return shapes

