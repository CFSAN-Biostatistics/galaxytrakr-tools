import pytest

from bio2srr import *


def test_element_tree_xpath():
    from xml.etree import ElementTree as xml
    root = xml.fromstring(biosample_example)
    assert root.find(".//Id[@db='BioSample']") is not None

def test_flatten_biosample_xml():
    d = flatten_biosample_xml(biosample_example)
    assert d['biosample_accession'] == 'SAMN17131268'
    assert d['organism'] == 'Campylobacter jejuni'
    assert d['isolate'] == 'CFSAN091032'

def test_flatten_runs():
    d = list(flatten_runs(runs_example))
    assert len(d) == 2

def test_header_sort_override_consistency():
    import random
    L = ["C", "B", "A", "taxid", "bioproject"]
    L.sort(key=hso)
    # assert L[0] == "bioproject"
    A = L.copy()
    assert A == L
    R = []
    for _ in range(100):
        random.shuffle(A)
        A.sort(key=hso)
        R.append(A == L)
    assert all(R)

def test_hso_override():
    assert header_sort_override("bioproject", "taxid") < 0
    assert header_sort_override("taxid", "bioproject") > 0
    assert header_sort_override("taxid", "taxid") == 0

def test_hso_regular():
    assert header_sort_override("A", "B") < 0
    assert header_sort_override("B", "A") > 0