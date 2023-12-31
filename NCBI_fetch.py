#!/usr/bin/env python3

import http.client
import json
import logging
import sys
import time
from pathlib import Path
from time import sleep
from urllib.error import HTTPError

import keyring
import toml
from Bio import Entrez

logging.basicConfig(
    filename=f"./logs/Entrez_info_{int(time.time())}.log",
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

email = keyring.get_password("Entrez", "Entrez_email")
api_key = keyring.get_password("Entrez", "Entrez_apikey")

if email is None:
    email = input("Entrez email: ")
    keyring.set_password("Entrez", "Entrez_email", email)
    logging.info("Set new email to keyring")

if api_key is None:
    api_key = input("Entrez api key: ")
    keyring.set_password("Entrez", "Entrez_apikey", api_key)
    logging.info("Set new api_key to keyring")


Entrez.email = email
Entrez.api_key = api_key


def query_builder(terms):
    terms_str = " OR ".join([f"{x}[Title]" for x in terms])
    query = (
        "("
        + terms_str
        + ") "
        + "AND phage[Title] "
        + "NOT hypothetical[Title] "
        + "NOT putative[Title] "
        + "NOT putitive[Title] "
        + "NOT probable[Title] "
        + "NOT possible[Title] "
        + "NOT unknown[Title] "
        + "AND 50:1000000[SLEN]"
    )

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
    out_file,
    batch_size=500,
    start_batch=0,
    cls=None,
    ret_mode="fasta",
    ret_type="gp",
):
    count = int(esearch_handler["Count"])

    for start in range(start_batch, count, batch_size):
        logging.info(
            f"{cls} - start: {start}, end: {start + batch_size}, total: {count}"
        )

        attempt = 0
        while attempt < 20:
            try:
                fetch_handle = Entrez.efetch(
                    db="protein",
                    retmode=ret_mode,
                    rettype=ret_type,
                    retstart=start,
                    retmax=batch_size,
                    webenv=esearch_handler["WebEnv"],
                    query_key=esearch_handler["QueryKey"],
                    idtype="acc",
                )
                attempt = 0
                data = fetch_handle.read()
                fetch_handle.close()
                with open(out_file, "ab" if type(data) == bytes else "a") as out:
                    out.write(data)
                logging.info(
                    f"\t\tCurrent file size: {'{:.2f}'.format(Path(out_file).stat().st_size/1024**2)} MB"
                )
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
                    f"{cls} - start: {start} | Received urllib HTTP error. Attempt {attempt}"
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

        if attempt >= 20:
            logging.error("Reached max number of attempts in a row without success")
            raise HTTPError


def main():
    config = toml.load("positive_terms.toml")
    class_labels = config["labels"]

    for cls, terms in class_labels.items():
        out_file = Path(f"./data/{cls}.xml")
        if out_file.is_file():
            logging.info(f"File {out_file} already exists! Halting program.")
            raise FileExistsError

        out_file.parent.mkdir(parents=True, exist_ok=True)

        query = query_builder(terms)
        esearch_handler = get_search(query)
        get_sequences(esearch_handler, out_file, cls=cls, ret_mode="xml")


if __name__ == "__main__":
    sys.exit(main())
