#!/usr/bin/env python3

# See https://www.ensembl.org/info/website/upload/gff.html
gtf_header = 'seqname source feature start end score strand frame attribute'.split()

def usable(record):
    """ Helper function to determine if we can use a record

    We are only interested in protein coding exons.

    """
    return (record['feature'] == 'exon' and
            record['attribute']['transcript_biotype'] =='protein_coding')


def parse_attribute(attribute):
    """ attribute - A semicolon-separated list of tag-value pairs """
    attributes = dict()
    for attr in attribute.split(';')[:-1]:
        key, *value = attr.strip().split(' ')
        attributes[key] = ' '.join(val.replace('"','') for val in value)
    return attributes


def gtf_to_json(gtf_handle, chrom_prefix=''):
    """ Parse a GTF file, and yield the records one by one """
    for line in gtf_handle:
        # Skip the headers
        if line.startswith('#'):
            continue
        spline = line.strip('\n').split('\t')
        record = {k:v for k, v in zip(gtf_header, spline)}
        record['attribute'] = parse_attribute(record['attribute'])
        record['seqname'] = f"{chrom_prefix}{record['seqname']}"
        yield record


def json_to_gtf(record):
    """ Return a record in GTF format """

    # First, we must package the attributes
    attr = ''.join(f' {key} "{value}";' for key, value in record['attribute'].items())
    # Strip the leading space
    record['attribute'] = attr.strip(' ')
    return '\t'.join(record[field] for field in gtf_header)

def gtf_to_csv(record, attribute_fields):
    """ Return a record in CSV format """
    gtf_fields = [record[field] for field in gtf_header if field != 'attribute']
    attribute_fields = [record['attribute'].get(field, '') for field in attribute_fields]
    return ','.join(gtf_fields + attribute_fields)
