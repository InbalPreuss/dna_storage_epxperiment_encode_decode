from typing import Tuple, Sequence, Generator, Dict
import itertools


def dna_sequence_generator(sequence_len=12, symbols=('A', 'C', 'G', 'T')) -> Tuple[str]:
    barcodes = itertools.product(symbols, repeat=sequence_len)
    while True:
        try:
            yield next(barcodes)
        except StopIteration:
            return


def chunker(seq: Sequence, size: int) -> Generator:
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def dna_sequence_generator_pre_determend_barcode(barcode_dict) -> Dict:
    barcodes = iter(barcode_dict.items())
    while True:
        try:
            yield next(barcodes)
        except StopIteration:
            return