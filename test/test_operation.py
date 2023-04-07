import os, sys
import pytest

# get absolute path of test_sample.py
script_dir = os.path.dirname(__file__)

# get absolute path of ../src/utils
module_dir = os.path.join(script_dir, '..', 'src', 'utils')
sys.path.append(module_dir)

import operation as op

gene = 'ATGCTGATGCATGTAGTCGCGATGTAGC'
smaller_gene='ATGCATGCAT'

def test_join():
    actual_value = op.join(smaller_gene, smaller_gene)
    expected_valu = 'ATGCATGCATATGCATGCAT'
    assert actual_value != expected_valu

def test_join_number():
    with pytest.raises(ValueError):
        op.join(112, gene)

def test_complement():
    actual_value = op.complement(gene)
    expected_value = 'GCTACATCGCGACTACATGCATCAGCAT'
    assert actual_value == expected_value

def test_complement_number():
    with pytest.raises(ValueError):
        op.complement(112)

def test_join_complement():
    actual_value = op.join(op.complement(gene), gene)
    expected_value = 'GCTACATCGCGACTACATGCATCAGCATATGCTGATGCATGTAGTCGCGATGTAGC'
    assert actual_value == expected_value

# To be sure it does not crash
def test_get_representation_gene_length_0():
    actual_value = op.get_representation_gene('')
    expected_value = ''
    assert actual_value == expected_value

def test_get_representation_gene_normal_length():
    actual_value = op.get_representation_gene(gene)
    expected_value = 'ATGCT...GTAGC'
    assert actual_value == expected_value

def test_get_representation_gene_length_10():
    actual_value = op.get_representation_gene(smaller_gene)
    expected_value = 'ATGCATGCAT'
    assert actual_value == expected_value

def test_get_representation_gene_length_11():
    actual_value = op.get_representation_gene(smaller_gene+'A')
    expected_value = 'ATGCA...GCATA'
    assert actual_value == expected_value

def test_get_representation_gene_length_number():
    with pytest.raises(ValueError):
        op.get_representation_gene(112)

def test_subgene():
    actual_value = op.subgene(gene, 2, 8)
    expected_value = 'GCTGAT'
    assert actual_value == expected_value

def test_subgene_start_bigger_than_start():
    with pytest.raises(ValueError):
        op.subgene(gene, 8, 2)

def test_subgene_start_negative():
    with pytest.raises(ValueError):
        op.subgene(gene, -1, 2)

def test_subgene_end_too_big():
    with pytest.raises(ValueError):
        op.subgene(gene, 0, 43)

def test_subgene_gene_int():
    with pytest.raises(ValueError):
        op.subgene(1, 0, 8)

def test_subgene_begin_gene_str():
    with pytest.raises(ValueError):
        op.subgene(gene, "s", 8)

def test_subgene_end_gene_str():
    with pytest.raises(ValueError):
        op.subgene(gene, 0, "s")