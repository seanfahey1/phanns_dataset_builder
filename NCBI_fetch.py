#!/usr/bin/env python3

import logging
import re
import sys
import threading
from pathlib import Path
from time import sleep
from urllib.error import HTTPError

import toml
from Bio import Entrez

from logger import setup_logging
from utils import fetch_credentials, get_args, query_builder


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
        out_dir = args.data_dir / f"{cls}"
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
