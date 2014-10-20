import os, sys

import gzip

import pysam

from collections import namedtuple

def build_header(chrm_sizes_urls = [
        "http://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes",]):
    # if we are passed a single url, then make it into a list
    if isinstance(chrm_sizes_urls, str):
        chrm_sizes_urls = [chrm_sizes_urls,]
    
    import urllib.request, cgi
    chrm_lengths = {}
    for url in chrm_sizes_urls:
        # doesn't support with syntax...
        req = urllib.request.urlopen(url)
        _, params = cgi.parse_header(req.headers.get('Content-Type', ''))
        encoding = params.get('charset', 'utf-8')

        # build a dictionary of chromosome lengths
        url_chrm_lengths = ( line.decode(encoding).split()
                             for line in req
                             if line.strip() != "")
        for contig_name, length in url_chrm_lengths:
            if contig_name in chrm_lengths:
                assert int(length) == chrm_lengths[contig_name]
            chrm_lengths[contig_name] = int(length)
        
        req.close()

    header = { 'HD': {'VN': '1.0'},
               'SQ': [{'LN': length, 'SN': contig_name}
                      for contig_name, length in chrm_lengths.items()] 
               }
    
    return header

class TagAlignRead(namedtuple('tagAlignLine', [
                   'contig', 'start', 'stop', 'seq', 'score', 'strand'])):
    pass

class TagAlign:
    def __init__(self, fname):
        if fname.endswith(".gz"):
            self.fp = gzip.open(fname, mode='rt')
        else:
            self.fp = open(fname)
        return
    
    def __enter__(self):
        return self
    
    def __iter__(self):
        self.fp.seek(0)
        for line in self.fp:
            data = line.split()
            data[1] = int(data[1])
            data[2] = int(data[2])
            data[4] = int(data[4])
            assert data[5] in "+-."
            yield TagAlignRead(*data)
        return
        
    def __exit__(self, type, value, tb):
        self.fp.close()

class OutputSam(pysam.Samfile):
    def add_tagalign_read(self, ta_rd, qname):
        rd_len = ta_rd.stop - ta_rd.start
        rd = pysam.AlignedRead()
        rd.qname = qname
        rd.seq = ( ta_rd.seq 
                   if len(ta_rd.seq) == rd_len
                   else "N"*rd_len )
        rd.flag = 1 + (16 if ta_rd.strand == '-' else 0)
        rd.rname = self.gettid(ta_rd.contig)
        rd.pos = ta_rd.start
        rd.mapq = ta_rd.score
        rd.cigar = ( (0, rd_len), )
        rd.mpos = -1
        rd.mrnm = -1
        #rd.tags = ( ("NM", 1),
        #           ("RG", "L1") )
        self.write(rd)

def main():
    with OutputSam("-", "wh", header=build_header()) as ofp:
        with TagAlign(sys.argv[1]) as inputfile:
            for i, rd in enumerate(inputfile):
                ofp.add_tagalign_read(rd, str(i))

if __name__ == '__main__':
    main()
