import contextlib

import click
import joblib
from joblib import Parallel, delayed
from tqdm import tqdm

import scphylo as scp


@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """[summary].

    Parameters
    ----------
    tqdm_object : [type]
        [description]

    Returns
    -------
    [type]
        [description]

    Yields
    -------
    [type]
        [description]
    """

    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()


def function(patient, tqdm_pos, kind):
    """[summary].

    Parameters
    ----------
    patient : [type]
        [description]
    tqdm_pos : [type]
        [description]
    """
    # for _ in tqdm(
    #     range(10), ascii=True, ncols=60, desc=f"", leave=False, position=pos + 1
    # ):
    #     sleep(1)
    scp.io.build_crc(patient=patient, tqdm_pos=tqdm_pos, kind=kind)


@click.command(short_help="Building the matrices of CRCs.")
def build_crcs():
    """[summary].

    Returns
    -------
    [type]
        [description]
    """
    patients = ["CRC01", "CRC02", "CRC04", "CRC09", "CRC10", "CRC11"]
    kinds = ["FPKM", "FPKM", "TPM", "TPM", "TPM", "TPM"]

    # only met
    # "CRC12", "CRC13", "CRC14", "CRC15"

    number_of_samples = len(patients)
    number_of_threads = len(patients)
    with tqdm_joblib(
        tqdm(
            ascii=True,
            ncols=75,
            desc="STATUS",
            total=number_of_samples,
            position=0,
        )
    ) as _:
        Parallel(n_jobs=number_of_threads)(
            delayed(function)(patients[i], i + 1, kinds[i])
            for i in range(number_of_samples)
        )

    # function("CRC01", 1, "FPKM")

    return None
