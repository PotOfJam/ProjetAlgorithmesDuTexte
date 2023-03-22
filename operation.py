def join(first, second):
    """Join operator for the eukaroytes

    Args:
        first (int): First part to join
        second (int): Second part to join

    Returns:
        int: The combination of both parts
    """
    return str(first) + str(second)


def complement(gene):
    """Complement operator for the eukaryotes

    Args:
        gene (string): The gene for which we want to obtain the complement

    Returns:
        string: The complement of the gene
    """
    res = ""
    for i in range(len(gene)):
        changed = False
        if (gene[i] == 'A'):
            res = res + 'T'
            changed = True
        elif (changed == False and gene[i] == 'T'):
            res = res + 'A'
            changed = True
        elif (changed == False and gene[i] == 'G'):
            res = res + 'C'
            changed = True
        elif (changed == False and gene[i] == 'C'):
            res = res + 'G'
            changed = True
    return res


def get_representation_gene(gene):
    """Visual representation of a gene:
    five first characters of the gene...five last characters of the gene

    Args:
        gene (string): The gene we want to obtain the representation

    Returns:
        string: The gene itself or its representation
    """
    return gene[0:5] + "..." + gene[-5:] if (len(gene) > 10) else gene


def subgene(gene, begin_gene, end_gene):
    """Extract part of a gene

    Args:
        gene (string): The gene for which we want to extract a subpart
        begin_gene (int): The begin index of the subpart we want to extract
        end_gene (int): The end index of the subpart we want to extract

    Raises:
        ValueError: When the begin_gene or end_gene does not respect the if-conditions

    Returns:
        string: The subpart of the gene
    """
    gene_length = len(gene)
    if (begin_gene < gene_length-1 and begin_gene >= 0 and \
                                        end_gene <= gene_length-1 and \
                                        end_gene > 0 and \
                                        begin_gene < end_gene):
        return gene[begin_gene:end_gene];
    else :
        raise ValueError("Error: begin_gene / end_gene incorrect.")



# gene = 'ATGCTGATGCATGTAGTCGCGATGTAGC'
# smaller_gene='ATGCATGCAT'
# print(join(complement(gene), gene))
# print(get_representation_gene(gene))
# print(get_representation_gene(smaller_gene))
# print(get_representation_gene(smaller_gene+'A'))
# print(subgene(gene, 2, 8))
# # print(subgene(gene, 8, 2)) # Should raise a ValueError!
# # print(subgene(gene, -1, 2)) # Should raise a ValueError!
# # print(subgene(gene, 0, 43)) # Should raise a ValueError!