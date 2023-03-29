def join(first, second):
    """Join operator for the eukaroytes

    Args:
        first (string): First part to join
        second (string): Second part to join

    Returns:
        int: The combination of both parts
    """
    if type(first) is not str or type(second) is not str:
        raise ValueError("Error: one of the gene is not a string.")

    return first + second


def complement(gene):
    """Complement operator for the eukaryotes

    Args:
        gene (string): The gene for which we want to obtain the complement

    Returns:
        string: The complement of the gene
    """

    if type(gene) is not str:
        raise ValueError("Error: one of the gene is not a string.")

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
    return res[::-1]


def get_representation_gene(gene):
    """Visual representation of a gene:
    five first characters of the gene...five last characters of the gene

    Args:
        gene (string): The gene we want to obtain the representation

    Returns:
        string: The gene itself or its representation
    """
    if type(gene) is not str:
        raise ValueError("Error: gene should be a string.")

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

    if type(gene) is not str or type(begin_gene) is not int or type(end_gene) is not int:
        raise ValueError("Error: one of the parameter has an incorrect type.")

    gene_length = len(gene)
    if (begin_gene < gene_length-1 and begin_gene >= 0 and \
                                        end_gene <= gene_length-1 and \
                                        end_gene > 0 and \
                                        begin_gene < end_gene):
        return gene[begin_gene:end_gene]
    else :
        raise ValueError("Error: begin_gene / end_gene incorrect.")