"""
Tests for `vcffuse` module.
"""
import pytest
from vcffuse import CreateSVGPicture, ExonCoords
from intervaltree import Interval, IntervalTree

class TestVcffuse(object):

    @classmethod
    def setup_class(cls):
        pass

    def test_something(self):
        pass

    def test_SVG_creation(self):
        exon_intervals = ( IntervalTree(
            [Interval(0, 300),
             Interval(400, 500),
             Interval(600, 700)]),
         IntervalTree(
            [Interval(800, 900),
             Interval(1000, 1100)]))
        SVG_creator = CreateSVGPicture.CreateSVGPicture()
        SVG_creator.verbose = True
        assert True == SVG_creator.verbose
        counter = SVG_creator.make_SVG_file(exon_intervals, "SVG_creator_test.svg", 9, "PRIME5", "PRIME3")
        assert 10 == counter

    @classmethod
    def teardown_class(cls):
        pass
