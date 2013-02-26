#!/usr/bin/env python
"""Script to inject into QIIME (1.3) workflow that removes discordant
read-pairs from OTU maps produced in two separate runs (one for each
end). OTUs must have same ids, i.e. this only works if you used a OTU
reference. Output is a QIIME-style OTU table."""


#--- standard library imports
#
import logging
import os
import sys
from optparse import OptionParser
import gzip

#--- third-party imports
#
# invocation of ipython on exceptions
if False:
    import sys, pdb
    from IPython.core import ultratb
    sys.excepthook = ultratb.FormattedTB(mode='Verbose',
                                         color_scheme='Linux', call_pdb=1)

#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = ""
__license__ = ""
__credits__ = [""]
__status__ = ""


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')



def build_read_id_map(fname_fa):
    """
    FIXME:add-doc
    """
    fh_fa = open(fname_fa, 'r')
    read_id_map = dict()
    for line in fh_fa:
        line = line.rstrip(os.linesep)
        if len(line) == 0 or line.startswith("#"):
            continue

        if not line.startswith(">"):
            continue
        line_split = line[1:].split()
        if len(line_split)!=2:
            LOG.info("Hope you know what you're doing. Was expecting exactly two elements on read line '%s'" % line)
        (qiime_id, orig_id) = line_split[:2]
        assert read_id_map.has_key(qiime_id) == False, (
            "Encountered non-unique QIIME read-id '%s'" % qiime_id)
        read_id_map[qiime_id] = orig_id
    fh_fa.close()
    
    LOG.info("Read %d ids from %s" % (len(read_id_map), fname_fa))
    return read_id_map



def parse_otu_tax(fname_tax):
    """
    Needs to be greengenes like input (easier parsing), i.e.
    OTU id\tk__Bacteria;p__...;c__...;o__...;f__...;g__...;s_...
    """

    otu_tax_map = dict()
    if fname_tax[-3:] == ".gz":
        fh_tax = gzip.open(fname_tax, 'r')
    else:
        fh_tax = open(fname_tax, 'r')
    for line in fh_tax:
        line = line.rstrip(os.linesep)
        if len(line) == 0 or line.startswith("#"):
            continue

        line_split = line.split('\t')
        assert len(line_split)==2, (
            "Expected exactly two values on the following line"
            " from '%s': '%s'" % (
                line, fname_tax))
        (otu, tax) = line_split
        otu_tax_map[otu] = tax
        assert tax.startswith("k__"), (
            "Non-Greengenes beginning of taxon '%s' in file '%s'" % (
                tax, fname_tax))

    LOG.info("Read %d OTU classifications from %s" % (
        len(otu_tax_map), fname_tax))
    return otu_tax_map
    

    
def parse_otu_map(fname_otu, read_id_map):
    """Use value of read_id_map[qiime_id] instead of qiime_id as key here
    and otu as value
    """

    IGNORE_NONE = False
    # if you ignore none's the number of fp goes up
    
    qiime_read_id_to_otu_map = dict()
    fh_otu = open(fname_otu, 'r')
    for line in fh_otu:
        line = line.rstrip(os.linesep)
        if len(line) == 0 or line.startswith("#"):
            continue

        line_split = line.split()
        otu = line_split[0]
        if IGNORE_NONE and otu.startswith('None'):
            LOG.debug("Ignoring assignments to"
                     " unknown OTU in line '%s'" % (line))
            continue
        for qiime_id in line_split[1:]:
            try:
                orig_id = read_id_map[qiime_id]
            except KeyError:
                sys.stderr.write("WARNING: orig_id = read_id_map[qiime_id] failed. Dropping to console.")
                import pdb; 
                pdb.set_trace()
            qiime_read_id_to_otu_map[orig_id] = otu
    fh_otu.close()
    
    LOG.info("Read assignments to %d OTUs from %s" % (
        len(set(qiime_read_id_to_otu_map.values())), fname_otu))
    return qiime_read_id_to_otu_map



def filter_and_join_paired_reads(method,
                                 f_seqs_1, f_otu_map_1,
                                 f_seqs_2, f_otu_map_2,
                                 f_tax, f_out):
    """FIXME:add-doc
    """
    
    # qiime id's are different to the original ones and the ones
    # referring to the same read pair might be different, since qiime
    # does filtering etc and then renames the ids according to sample
    # name and a running number

    read_id_map_1 = build_read_id_map(f_seqs_1)
    read_id_map_2 = build_read_id_map(f_seqs_2)
    #import pdb; pdb.set_trace()
    otu_map_1 = parse_otu_map(f_otu_map_1, read_id_map_1) 
    otu_map_2 = parse_otu_map(f_otu_map_2, read_id_map_2)
    
    otu_tax_map = parse_otu_tax(f_tax)

    # loop over id's present in both otu mappings
    same_otu = 0
    diff_otu = 0

    otu_map_out = dict()
    
    shared_orig_ids = set(otu_map_1.keys()).intersection(set(otu_map_2.keys()))
    LOG.info("Will check %d (original) read-pair ids shared between both" % len(shared_orig_ids))
    for orig_id in shared_orig_ids:
        otu_id_1 = otu_map_1[orig_id]
        otu_id_2 = otu_map_2[orig_id]

        if otu_id_1 == otu_id_2:
            same_otu += 1
        else:
            diff_otu += 1
            
        if otu_id_1.startswith("None"):
            tax_otu_1 = 'None'
        else:
            tax_otu_1 = otu_tax_map[otu_id_1]

        if otu_id_2.startswith("None"):
            tax_otu_2 = 'None'
        else:
            tax_otu_2 = otu_tax_map[otu_id_2]

            
        # method: otu must be the same between pairs?
        #
        if method == 'same-otu':
            if otu_id_1 == otu_id_2:
                # use 1 or 2...they are the same
                otu_map_out[otu_id_1] = otu_map_out.get(otu_id_1, 0) + 1
            else:
                LOG.debug("discard %s: %s vs %s" % (
                    orig_id, tax_otu_1, tax_otu_2))
                continue

        # method: otu must have same taxonomy?
        #
        elif method == 'same-tax':
            if tax_otu_1 == tax_otu_2:
                # keep the otu which already has the most and the
                # first one if equal
                if otu_map_out.get(otu_id_1, 0) >= otu_map_out.get(otu_id_2, 0):
                    otu_map_out[otu_id_1] = otu_map_out.get(otu_id_1, 0) + 1
                else:
                    otu_map_out[otu_id_2] = otu_map_out.get(otu_id_2, 0) + 1
                    
            else:
                LOG.debug("discard %s: %s vs %s" % (
                    orig_id, tax_otu_1, tax_otu_2))
                continue

            
        # method: concordant taxonomy, i.e. they don't contradict
        #
        elif method == 'concordant-tax':
            # comparison of taxonomic assignment. if they don't
            # contradict, keep most detailed taxon.
            
            # remove everything up to first ambigious level Either
            # empty (__;) or starts with 'unclassified' and
            # 'uncultured' in Swee Hoe's DB
            for ambiguity_string in  ['__;', '__unc']:
                try:
                    tax_otu_1_cmp = tax_otu_1[:tax_otu_1.index(ambiguity_string)-1]
                except ValueError:
                    tax_otu_1_cmp = tax_otu_1
                try:
                    tax_otu_2_cmp = tax_otu_2[:tax_otu_2.index(ambiguity_string)-1]
                except ValueError:
                    tax_otu_2_cmp = tax_otu_2
     
                ambiguity_string = ';unc'
                try:
                    tax_otu_1_cmp = tax_otu_1[:tax_otu_1.index(ambiguity_string)]
                except ValueError:
                    tax_otu_1_cmp = tax_otu_1
                try:
                    tax_otu_2_cmp = tax_otu_2[:tax_otu_2.index(ambiguity_string)]
                except ValueError:
                    tax_otu_2_cmp = tax_otu_2


            
            min_len = min([len(tax_otu_1_cmp), len(tax_otu_2_cmp)])
            if tax_otu_1_cmp[:min_len] == tax_otu_2_cmp[:min_len]:
                LOG.debug("match: %s vs %s" % (tax_otu_1, tax_otu_2))
                if len(tax_otu_1) >= len(tax_otu_2):
                    most_spec_otu = otu_id_1
                else:
                    most_spec_otu = otu_id_2
                otu_map_out[most_spec_otu] = otu_map_out.get(most_spec_otu, 0) + 1
            else:
                LOG.debug("discard %s: %s vs %s" % (
                    orig_id, tax_otu_1, tax_otu_2))
                continue

        # method: unknown
        #
        else:
            raise ValueError, ("Unknown method")

    LOG.info("%d OTUs were identical between pairs; %d differed." % (
        same_otu, diff_otu))


    fh_out = open(f_out, 'w')
    fh_out.write("# QIIME v1.3.0 OTU table\n")
    fh_out.write("#OTU ID\tFIXME:samplename\tConsensus Lineage\n")
    for (otu_id, count) in sorted(otu_map_out.items(),
                                         key=lambda x: x[1]):
        if otu_id.startswith("None"):
            lineage = 'None'
        else:
            lineage = otu_tax_map[otu_id]
        fh_out.write("%s\t%d\t%s\n" % (otu_id, count, lineage))
    fh_out.close()

    LOG.info("Written OTU table '%s'. Please fix header if neccessary" % f_out)



def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog [options]\n\n" + __doc__
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="debugging")
    parser.add_option("-o", "--output",
                      dest="fout", # type="string|int|float"
                      help="output OTU table")
    parser.add_option("-t", "--taxonomy",
                      dest="ftax", # type="string|int|float"
                      help="Greengenes-style taxonomy assignment file classifying seqs/OTUs")
    parser.add_option("", "--seqs1",
                      dest="fseq1", # type="string|int|float"
                      help="Fasta file containing the QIIME formatted reads of one end"
                      ", i.e. the seq ids are the qiime-id followed by the original read id")
    parser.add_option("", "--seqs2",
                      dest="fseq2", # type="string|int|float"
                      help="Fasta file containing the QIIME formatted reads of the other end"
                      ", i.e. the seq ids are the qiime-id followed by the original read id")
    parser.add_option("", "--otu1",
                      dest="fotu1", # type="string|int|float"
                      help="OTU map (i.e., result from pick_otus.py, NOT otu_table.txt) for one end")
    parser.add_option("", "--otu2",
                      dest="fotu2", # type="string|int|float"
                      help="OTU map (i.e., result from pick_otus.py, NOT otu_table.txt) for other end")
    choices = ['same-otu', 'same-tax', 'concordant-tax']
    default = "same-otu"
    parser.add_option("-m", "--method",
                      dest="method",  
                      default=default,
                      choices=choices,
                      help="OTU map (i.e., result from pick_otus.py) for other end")

    return parser



def main():
    """FIXME:add-doc
    """

    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    if len(args):
        parser.error("Unrecognized arguments found: %s." % (
            ' '.join(args)))
        sys.exit(1)
        
    if opts.verbose:
        LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)
        
    if not opts.method:
        parser.error("Missing argument: method.")
        sys.exit(1)
    if not opts.fotu1 or not opts.fotu2:
        parser.error("Missing argument: One or both OTU map files.")
        sys.exit(1)
    if not opts.fseq1 or not opts.fseq2:
        parser.error("Missing argument: One or both fasta sequence files.")
        sys.exit(1)
    if not opts.ftax:
        parser.error("Missing argument: Taxonomy assignment file.")
        sys.exit(1)
    
    for file_to_check in [opts.fotu1, opts.fotu2,
                          opts.fseq1, opts.fseq2, opts.ftax]:
        if not os.path.exists(file_to_check):
            LOG.critical("Non-existant file '%s'." % file_to_check)
            sys.exit(1)

    if not opts.fout:
        parser.error("Missing argument: output OTU table.")
        sys.exit(1)
    if os.path.exists(opts.fout):
        LOG.critical("Cowardly refusing to"
                     " overwrite existant file '%s'." % opts.fout)
        sys.exit(1)

        
    if opts.method == 'same-otu':
        LOG.warn("If same-otu becomes the default, then this script can be simplified (no tax needed, otu map instead of table output etc")

    filter_and_join_paired_reads(opts.method, opts.fseq1, opts.fotu1,
                                 opts.fseq2, opts.fotu2,
                                 opts.ftax, opts.fout)

          

if __name__ == "__main__":
    main()
    LOG.info("successful exit")
