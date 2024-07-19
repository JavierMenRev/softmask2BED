import argparse
from Bio import SeqIO
from sys import stdin as std_in, stdout as std_out

def find_lowercase_stretches(seq):
    stretches = []
    in_stretch = False
    start = None
    
    for i, nucleotide in enumerate(seq):
        if nucleotide.islower():
            if not in_stretch:
                start = i
                in_stretch = True
        else:
            if in_stretch:
                end = i - 1
                stretches.append((start, end))
                in_stretch = False
    
    if in_stretch:
        stretches.append((start, len(seq) - 1))
    
    return stretches

def get_options():
    parser = argparse.ArgumentParser(description='Convert soft-masked regions in a FASTA file to BED format.')
    parser.add_argument('-I', metavar='input', nargs='?', default=std_in, required=False,
                        help='Path to the input FASTA file with soft-masked nucleotides, defaults to standard input.')
    parser.add_argument('-O', metavar='output', nargs='?', default=std_out, required=False,
                        help='Path to the output BED file, defaults to standard output.')
    return parser.parse_args()

def main():
    args = get_options()
    
    input_file = args.I if args.I != std_in else None
    output_file = args.O if args.O != std_out else None

    with (open(output_file, "w") if output_file else std_out) as bed_file:
        for record in SeqIO.parse(input_file, "fasta"):
            seq = str(record.seq)
            stretches = find_lowercase_stretches(seq)
            for start, end in stretches:
                bed_file.write(f"{record.id}\t{start}\t{end + 1}\n")

if __name__ == '__main__':
    main()
