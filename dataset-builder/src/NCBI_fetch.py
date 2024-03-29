#!/usr/bin/env python3

import argparse
import logging
import re
import sys
import threading
import time
from pathlib import Path
from time import sleep
from urllib.error import HTTPError

import keyring
import toml
from Bio import Entrez


def setup_logging(args, directory):
    log_dir = Path(directory)
    log_dir.mkdir(exist_ok=True, parents=True)

    logging.basicConfig(
        filename=log_dir / f"Entrez_info_{int(time.time())}.log",
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    for arg, value in args.__dict__.items():
        logging.info(f"{arg}: {value}")


def get_args():
    parser = argparse.ArgumentParser(
        prog="NCBI fetch",
        description="""
Query the NCBI protein database using Entrez
This tool requires a configuration toml file in the following format:

[positive_labels]
LABEL1 = ["term1", "term2", ... "termN"]
LABEL2 = ["term1", "term2", ... "termN"]
...
LABELN = ["term1", "term2", ... "termN"]

[query]
additional_query = [
    "AND refseq[filter]",
    "AND phage[Title]",
    "NOT hypothetical[Title]",
    "NOT putative[Title]",
    "NOT putitive[Title]",
    "NOT probable[Title]",
    "NOT possible[Title]",
    "NOT unknown[Title]",
    "AND 50:1000000[SLEN]"
    ]

""",
        epilog="Example command: \npython NCBI_fetch.py -c path/to/config/file.toml -e [you@abc.edu] -a [api key]",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-c", "--config", required=True, help="Path to (the configuration toml file."
    )
    parser.add_argument(
        "-e", "--email", nargs="?", help="Entrez email address (username)."
    )
    parser.add_argument("-a", "--api_key", nargs="?", help="Entrez API key.")
    parser.add_argument(
        "-d",
        "--data_dir",
        default="./data",
        help="Relative path to directory to store data (will be created if it doesn't exist).",
    )
    parser.add_argument(
        "-l",
        "--log_dir",
        default="./logs",
        help="Relative path to directory to store log file (will be created if it doesn't exist).",
    )

    args = parser.parse_args()
    return args


def fetch_credentials(email, api_key):
    if email is None:
        email = keyring.get_password("Entrez", "Entrez_email")
        if email is None:
            email = input("Entrez email: ")
            keyring.set_password("Entrez", "Entrez_email", email)
            logging.info("Set new email to keyring")

    if api_key is None:
        api_key = keyring.get_password("Entrez", "Entrez_apikey")
        if api_key is None:
            api_key = input("Entrez api key: ")
            keyring.set_password("Entrez", "Entrez_apikey", api_key)
            logging.info("Set new api_key to keyring")

    Entrez.email = email
    Entrez.api_key = api_key
    return email, api_key


def query_builder(terms, additional_query):
    terms_str = " OR ".join([f"{x}[Title]" for x in terms])
    query = "(" + terms_str + ") " + " ".join(additional_query)

    logging.info(f"Built query as: {query}")

    return query


def get_search(query):
    with Entrez.esearch(
        db="protein", term=query, idtype="acc", usehistory="y"
    ) as handle:
        esearch_handler = Entrez.read(handle)

    logging.info(f'ESearch returned {esearch_handler.get("Count")} results')
    return esearch_handler


def get_sequences(
    esearch_handler,
    out_dir,
    batch_size=1,
    start_batch=0,
    cls=None,
    ret_mode="text",
    ret_type="fasta",
):
    count = int(esearch_handler["Count"])

    for start in range(start_batch, count, batch_size):
        logging.info(
            f"{cls} - start: {start}, end: {start + batch_size}, total: {count}"
        )

        attempt = 0
        while attempt < 100:
            try:
                with Entrez.efetch(
                    db="protein",
                    retmode=ret_mode,
                    rettype=ret_type,
                    retstart=start,
                    retmax=batch_size,
                    webenv=esearch_handler["WebEnv"],
                    query_key=esearch_handler["QueryKey"],
                    idtype="acc",
                ) as fetch_handle:
                    data = fetch_handle.read()

                name_match = re.search(">(.*?) ", data).group(1)
                base_name = name_match + "." + str(ret_type)
                out_file = out_dir / base_name
                if out_file.is_file():
                    logging.info(
                        f"File already exists with name {out_file}... Skipping."
                    )

                else:
                    with open(out_file, "wb" if type(data) == bytes else "w") as out:
                        out.write(data)
                logging.info(
                    f"\t\tCurrent number of files: {len(list(Path(out_dir).glob('*')))}"
                )
                attempt = 0
                break

            except HTTPError as err:
                attempt += 1
                logging.error(
                    f"{cls} - start: {start} | Received HTTP error. Attempt number {attempt}"
                )
                logging.error(err)
                sleep(15 * attempt)

            except ValueError as err:
                attempt += 1
                logging.error(
                    f"{cls} - start: {start} | Received urllib HTTP error or ValueError. Attempt {attempt}"
                )
                logging.error(err)
                sleep(180)

            except Exception as err:
                attempt += 1
                logging.error(
                    f"UNCAUGHT EXCEPTION | {cls} - start: {start} | Attempt number {attempt}"
                )
                logging.error(err)
                sleep(15 * attempt)

        if attempt >= 50:
            logging.error("Reached max number of attempts in a row without success")
            raise HTTPError


def main():
    args = get_args()
    setup_logging(args, args.log_dir)
    fetch_credentials(args.email, args.api_key)

    config = toml.load("config.toml")
    class_labels = config["positive_labels"]

    job_queue = []
    for cls, terms in class_labels.items():
        out_dir = Path(args.data_dir) / f"{cls}"
        out_dir.mkdir(parents=True, exist_ok=True)

        query = query_builder(terms, config["query"].get("additional_query"))
        esearch_handler = get_search(query)
        job_queue.append(
            threading.Thread(
                target=get_sequences,
                args=(
                    esearch_handler,
                    out_dir,
                    1,
                    0,
                    cls,
                    "text",
                    "fasta",
                ),
            )
        )

    for thread in job_queue:
        thread.start()

    for thread in job_queue:
        thread.join()


if __name__ == "__main__":
    sys.exit(main())
