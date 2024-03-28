import logging
import time
from pathlib import Path


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
