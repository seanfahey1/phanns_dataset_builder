import argparse
import logging

import keyring
from Bio import Entrez


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
