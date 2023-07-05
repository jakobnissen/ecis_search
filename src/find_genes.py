import Bio.SeqIO
import pyrodigal

#import os
#os.chdir("..")

records = list(Bio.SeqIO.parse("raw/toxin_hosts.fna", "fasta"))

with open("/tmp/foo.fna", "w") as file:
    for record in records:
        orf_finder = pyrodigal.OrfFinder()
        orf_finder.train(bytes(record.seq))
        genes = orf_finder.find_genes(bytes(record.seq))
        genes.write_genes(file, sequence_id=record.id, width=80)