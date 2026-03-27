from bz2 import open as bzopen
from gzip import open as gzopen

from contextlib import ExitStack
from itertools import zip_longest
from pathlib import Path
from sys import argv

import random


usage = """

"""

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

# file compression signatures
magics = {
    b'\x1f\x8b\x08':gzopen,
    b'\x42\x5a\x68':bzopen,
}

def sniff(path):
    "Sniff first three bytes of the file to determine format based on the magic number."
    with open(path, 'rb') as fp:
        magic = fp.read(3)
    return magics.get(magic, open)


def coverage(collection, genome_size):
    "Collection of 1 or 2 tuples, whose 2nd item is the read string"
    return sum((len(read[0][1]) for read in collection)) / genome_size # reverse read pair doesn't contribute to coverage so we can ignore it


try:
    fin, rin, fout, rout, cov, gen_size, *opts = argv[1:]
    ins = [fin, rin]
    outs = [fout, rout]
except ValueError: # not enough values to unpack
    try:
        fin, fout, cov, gen_size, *opts = argv[1:]
        ins = [fin]
        outs = [fout]
    except ValueError:
        print(usage)
        quit(1)
try:
    cov = float(cov)
    gen_size = int(gen_size)
except ValueError:
    print("Desired coverage and assumed genome size should be numbers")
    print(usage)
    quit(1)

seed = "ed2b99d842cddc1ac81d7c01a0bf0555"
if opts:
    seed = opts[0]
random.seed(seed)

assert len(ins) == len(outs)
file_openers = [sniff(path) for path in ins] # output format determined by input format
with ExitStack() as stack:
    ins = [stack.enter_context(openn(path, 'r')) for openn, path in zip(file_openers, ins)] # opened input files
    inns = [iter(grouper(inn, 4)) for inn in ins] # stateful 4-ply iterator over lines in the input
    outs = [stack.enter_context(openn(path, 'w')) for openn, path in zip(file_openers, outs)] # opened output files

    for file in ins:
        if hasattr(file, "name"):
            print(file.name)

    # https://en.m.wikipedia.org/wiki/Reservoir_sampling

    reservoir = []
    # this is going to be 1 or 2-tuples of 4-tuples representing the 4 lines of the fastq file
    # we determine its current coverage (and thus its reservoir size) to fill it, which consumes reads
    # from the open files
    reads = 0
    for i, readpair in enumerate(zip(*inns)):
        reads += len(readpair[0][1])
        reservoir.append(readpair)
        if reads / gen_size > cov:
            break

    k = len(reservoir) # this is about how big the reservoir needs to be to get cov coverage
    #W = exp(log(random.random()) / k)

    random.shuffle(reservoir)

    print(f"{k} reads selected to achieve {coverage(reservoir, gen_size):.3f}X coverage.")

    # if the number of reads is too few to meet the coverage cutoff, then the iterators
    # should be exhausted and this won't run
    # this is essentially Algorithm L, as I understand it
    for i, readpair in enumerate(zip(*inns)):
        r = random.randint(0, i)
        if r < k:
            reservoir[r] = readpair

    for readpair in reservoir: # output the sampled reads
        for read, file in zip(readpair, outs):
            defline, read, spacer, quals = read
            file.write(defline)
            file.write(read)
            file.write(spacer)
            file.write(quals)

# [fp.close() for fp in ins]
# [fp.close() for fp in outs]