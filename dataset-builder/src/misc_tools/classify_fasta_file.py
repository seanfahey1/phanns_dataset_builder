#! /usr/env/python

import argparse
import sys
from collections import Counter
from pathlib import Path

import toml
from Bio import SeqIO
from tqdm import tqdm


def verify_config_file(file_path):
    file = Path(file_path)

    assert file.is_file(), f"File {file_path} does not exist"
    assert file.suffix == ".toml", f"File {file_path} is not a toml file"

    return file


def verify_fasta_file(file_path):
    file = Path(file_path)

    assert file.is_file(), f"File {file_path} does not exist"
    assert file.suffix in (".fasta", ".faa"), f"File {file_path} is not a fasta file"

    return file


def get_args():
    parser = argparse.ArgumentParser(description="Classify fasta file")
    parser.add_argument(
        "-i", "--input", help="Input fasta file", required=True, type=verify_fasta_file
    )
    parser.add_argument("-o", "--output", help="Output file", required=False)
    parser.add_argument(
        "-c",
        "--config",
        help="Config file",
        default="config.toml",
        type=verify_config_file,
    )

    return parser.parse_args()


def header_checker(header, config):
    if "hypothetical" in header.lower():
        return "hypothetical protein"

    else:
        tag = None
        for term, values in config.items():
            for value in values:
                if value.lower() in header.lower():
                    if tag is None:
                        tag = term
                    else:
                        tag = "multiple matches"
                        return tag

        if tag is None:
            return "OTH"
        else:
            return tag


def main():
    args = get_args()

    config = toml.load(args.config).get("positive_labels", {})
    if args.output is None:
        args.output = args.input.with_suffix(".SF_annotated.fasta")

    output = []
    with open(args.input, "r") as input_file:
        fasta_sequences = list(SeqIO.parse(input_file, "fasta"))
        tag_counter = Counter()

        for fasta in tqdm(fasta_sequences):
            header = fasta.description

            tag = header_checker(header, config)
            tag_counter[tag] += 1

            fasta.description = f"{header} --SF_annotation: {tag}"
            output.append(fasta)

    for k, v in tag_counter.items():
        print(f"{k}: {v:,}")

    SeqIO.write(output, args.output, "fasta")
    print(f"Wrote file to {args.output}")


if __name__ == "__main__":
    sys.exit(main())
