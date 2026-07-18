import contextlib
import datetime
import functools
import importlib.resources
import multiprocessing
import os
import shutil
import tempfile
import time

import joblib
import numpy as np

import scphylo as scp


def log_input(df_in):
    """Log dimensions and genotype-state counts for an input matrix."""
    size = df_in.shape[0] * df_in.shape[1]
    scp.logg.info(f"input -- size: {df_in.shape[0]}x{df_in.shape[1]}")
    scp.logg.info(
        f"input -- 0: {np.sum(df_in.values == 0)}#,"
        f" {100 * np.sum(df_in.values == 0) / size:.1f}%"
    )
    scp.logg.info(
        f"input -- 1: {np.sum(df_in.values == 1)}#,"
        f" {100 * np.sum(df_in.values == 1) / size:.1f}%"
    )
    scp.logg.info(
        f"input -- NA: {np.sum(df_in.values == 3)}#,"
        f" {100 * np.sum(df_in.values == 3) / size:.1f}%"
    )
    scp.logg.info(f"input -- CF: {is_conflict_free_gusfield(df_in)}")


def log_output(df_out, running_time):
    """Log dimensions, genotype states, validity, and runtime for an output."""
    size = df_out.shape[0] * df_out.shape[1]
    scp.logg.info(f"output -- size: {df_out.shape[0]}x{df_out.shape[1]}")
    scp.logg.info(
        f"output -- 0: {np.sum(df_out.values == 0)}#,"
        f" {100 * np.sum(df_out.values == 0) / size:.1f}%"
    )
    scp.logg.info(
        f"output -- 1: {np.sum(df_out.values == 1)}#,"
        f" {100 * np.sum(df_out.values == 1) / size:.1f}%"
    )
    scp.logg.info(
        f"output -- NA: {np.sum(df_out.values == 3)}#,"
        f" {100 * np.sum(df_out.values == 3) / size:.1f}%"
    )
    icf = is_conflict_free_gusfield(df_out)
    scp.logg.info("output -- CF: ", end="")
    if icf:
        scp.logg.info(icf, color="green")
    else:
        scp.logg.info(icf, color="red")
    scp.logg.info(
        f"output -- time: {running_time:.1f}s"
        f" ({datetime.timedelta(seconds=running_time)})"
    )


def log_flip(df_in, df_out):
    """Log genotype flips and inferred error rates between two matrices."""
    flips_0_1, flips_1_0, flips_na_0, flips_na_1 = count_flips(
        df_in.values, df_out.values, 3
    )
    fn_rate, fp_rate, na_rate = infer_rates(df_in.values, df_out.values, 3)
    scp.logg.info(f"flips -- #0->1: {flips_0_1}")
    scp.logg.info(f"flips -- #1->0: {flips_1_0}")
    scp.logg.info(f"flips -- #NA->0: {flips_na_0}")
    scp.logg.info(f"flips -- #NA->1: {flips_na_1}")
    scp.logg.info(f"rates -- FN: {fn_rate:.3f}")
    scp.logg.info(f"rates -- FP: {fp_rate:.8f}")
    scp.logg.info(f"rates -- NA: {na_rate:.3f}")


def calc_nll_matrix(df_in, df_out, alpha, beta):
    """Calculate the negative log-likelihood of an inferred genotype matrix."""
    if alpha == 0 or beta == 0:
        return None
    columns = np.intersect1d(df_in.columns, df_out.columns)
    indices = np.intersect1d(df_in.index, df_out.index)
    D = df_in.loc[indices, columns].values
    E = df_out.loc[indices, columns].values
    removedMutations = []
    objective = 0
    for j in range(len(columns)):
        numZeros = 0
        numOnes = 0
        for i in range(len(indices)):
            if D[i, j] == 0:
                numZeros += 1
                objective += np.log(beta / (1 - alpha)) * E[i, j]
            elif D[i, j] == 1:
                numOnes += 1
                objective += np.log((1 - beta) / alpha) * E[i, j]
        objective += numZeros * np.log(1 - alpha)
        objective += numOnes * np.log(alpha)
        if j in removedMutations:
            objective -= numZeros * np.log(1 - alpha) + numOnes * (
                np.log(alpha) + np.log((1 - beta) / alpha)
            )
    return -objective


def stat(df_in, df_out, alpha, beta, running_time):
    """Log a summary comparing input and inferred genotype matrices."""
    log_input(df_in)
    log_output(df_out, running_time)
    log_flip(df_in, df_out)
    nll = calc_nll_matrix(df_in, df_out, alpha, beta)
    scp.logg.info(f"score -- NLL: {nll}")


def parse_params_file(filename):
    """Parse simulation and error-rate parameters encoded in a filename."""

    def _parse_params_file_helper(param):
        try:
            value = basename.split(f"{param}_")[1]
            if "-" in value:
                value = value.split("-")[0]
            else:
                value = value.split(".")[0]
            return float(value) if "." in value else int(value)
        except IndexError:
            return None

    data = {}
    _, basename = dir_base(filename)
    for param in [
        "simNo",
        "s",
        "m",
        "h",
        "minVAF",
        "ISAV",
        "n",
        "fp",
        "fn",
        "na",
        "d",
        "l",
    ]:
        value = _parse_params_file_helper(param)
        if value is not None:
            data[param] = value
    return data


def parse_log_file(filename):
    """Parse runtime and error-rate statistics from a solver log file."""
    result = {}
    _, basename = dir_base(filename)
    result["tool"] = basename.split(".")[-1]
    with open(filename) as fin:
        for line in fin:
            line = line.strip()
            if "output -- time: " in line:
                result["running_time"] = float(
                    line.replace("output -- time: ", "").split()[0].replace("s", "")
                )
            if "output -- CF: " in line:
                result["is_cf"] = bool(line.replace("output -- CF: ", "").split()[-1])
            if "rates -- FN: " in line:
                result["fn_rate"] = float(line.replace("rates -- FN: ", "").split()[-1])
            if "rates -- FP: " in line:
                result["fp_rate"] = float(line.replace("rates -- FP: ", "").split()[-1])
    return result


def parse_score_file(filename):
    """Parse key-value scores from a solver score file."""
    result = {}
    _, basename = dir_base(filename)
    result["tool"] = basename.split(".")[-1]
    with open(filename) as fin:
        for line in fin:
            line = line.strip()
            result[line.split("=")[0]] = float(line.split("=")[1])
    return result


def count_flips(I_mtr, O_mtr, na_value=3):
    """Count state changes between input and output genotype matrices."""
    flips_0_1 = 0
    flips_1_0 = 0
    flips_na_0 = 0
    flips_na_1 = 0
    n, m = I_mtr.shape
    for i in range(n):
        for j in range(m):
            if I_mtr[i, j] == 0 and O_mtr[i, j] == 1:
                flips_0_1 += 1
            elif I_mtr[i, j] == 1 and O_mtr[i, j] == 0:
                flips_1_0 += 1
            elif I_mtr[i, j] == na_value and O_mtr[i, j] == 0:
                flips_na_0 += 1
            elif I_mtr[i, j] == na_value and O_mtr[i, j] == 1:
                flips_na_1 += 1
    return flips_0_1, flips_1_0, flips_na_0, flips_na_1


def infer_rates(I_mtr, O_mtr, na_value=3):
    """Infer false-negative, false-positive, and missing-data rates."""
    flips_0_1, flips_1_0, flips_na_0, flips_na_1 = count_flips(I_mtr, O_mtr, na_value)
    fn_rate = flips_0_1 / ((O_mtr == 1) & (I_mtr != na_value)).sum()
    fp_rate = flips_1_0 / ((O_mtr == 0) & (I_mtr != na_value)).sum()
    na_rate = (flips_na_1 + flips_na_0) / I_mtr.size
    return fn_rate, fp_rate, na_rate


def is_conflict_free(df_in):
    """Check conflict-free criteria via Gusfield algorithm.

    The order of this algorithm is :math:`O(nm^2)`
    where n is the number of cells and m is the number of mutations.

    Parameters
    ----------
    df_in : :class:`pandas.DataFrame`
        Input genotype matrix.

    Returns
    -------
    :obj:`bool`
        A Boolean checking if the input conflict-free or not.

    See Also
    --------
    :func:`scphylo.ul.is_conflict_free_gusfield`.
    """
    D = df_in.astype(int).values
    if not np.array_equal(np.unique(D), [0, 1]):
        return False
    conflict_free = True
    for p in range(D.shape[1]):
        for q in range(p + 1, D.shape[1]):
            oneone = False
            zeroone = False
            onezero = False
            for r in range(D.shape[0]):
                if D[r, p] == 1 and D[r, q] == 1:
                    oneone = True
                if D[r, p] == 0 and D[r, q] == 1:
                    zeroone = True
                if D[r, p] == 1 and D[r, q] == 0:
                    onezero = True
            if oneone and zeroone and onezero:
                conflict_free = False
                return conflict_free
    return conflict_free


def is_conflict_free_gusfield(df_in):
    """Check conflict-free criteria via Gusfield algorithm.

    This is an implementation of algorithm 1.1 in :cite:`Gusfield_1991`.

    The order of this algorithm is :math:`O(nm)`
    where n is the number of cells and m is the number of mutations.

    Parameters
    ----------
    df_in : :class:`pandas.DataFrame`
        Input genotype matrix.

    Returns
    -------
    :obj:`bool`
        A Boolean checking if the input conflict-free or not.

    Examples
    --------
    >>> sc = scp.datasets.test()
    >>> scp.ul.is_conflict_free_gusfield(sc)
    False

    See Also
    --------
    :func:`scphylo.ul.is_conflict_free`.
    """
    I_mtr = df_in.astype(int).values
    if not np.array_equal(np.unique(I_mtr), [0, 1]):
        return False

    def _sort_bin(a):
        b = np.transpose(a)
        b_view = np.ascontiguousarray(b).view(
            np.dtype((np.void, b.dtype.itemsize * b.shape[1]))
        )
        idx = np.argsort(b_view.ravel())[::-1]
        c = b[idx]
        return np.transpose(c), idx

    Ip = I_mtr.copy()
    O_mtr, _ = _sort_bin(Ip)
    Lij = np.zeros(O_mtr.shape, dtype=int)
    for i in range(O_mtr.shape[0]):
        maxK = 0
        for j in range(O_mtr.shape[1]):
            if O_mtr[i, j] == 1:
                Lij[i, j] = maxK
                maxK = j + 1
    Lj = np.amax(Lij, axis=0)
    for i in range(O_mtr.shape[0]):
        for j in range(O_mtr.shape[1]):
            if O_mtr[i, j] == 1:
                if Lij[i, j] != Lj[j]:
                    return False
    return True


def tmpdir(prefix="scphylo.", suffix=".scphylo", dirname="."):
    """Create a temporary directory and return its path."""
    return tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dirname)


def tmpfile(prefix="scphylo.", suffix=".scphylo", dirname="."):
    """Create a temporary file and return its path."""
    return tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=dirname)[1]


def tmpdirsys(prefix="scphylo.", suffix=".scphylo", dirname="."):
    """Create a managed temporary-directory object."""
    return tempfile.TemporaryDirectory(suffix=suffix, prefix=prefix)


def cleanup(dirname):
    """Remove a directory tree."""
    shutil.rmtree(dirname)


def remove(filename):
    """Remove a file when it exists."""
    if os.path.exists(filename):
        os.remove(filename)


def dir_base(infile):
    """Return the directory and extension-free basename of a path."""
    basename = os.path.splitext(os.path.basename(infile))[0]
    dirname = os.path.dirname(infile)
    return dirname, basename


def dirbase(infile):
    """Return a path with its final extension removed."""
    return os.path.splitext(infile)[0]


def mkdir(indir):
    """Create a directory when absent and return its path."""
    if not os.path.exists(indir):
        os.makedirs(indir)
    return indir


def executable(binary, appname):
    """Locate an executable in PATH or the configured tools directory."""
    from scphylo.ul._external import resolve_executable

    return resolve_executable(binary, appname)


def timeit(f):
    """Wrap a callable to log its execution time."""

    def wrap(*args, **kwargs):
        start_time = time.time()
        ret = f(*args, **kwargs)
        end_time = time.time()
        if end_time - start_time < 60:
            scp.logg.info(f"Time needed for {f.__name__}: {end_time - start_time:.3f}")
        else:
            scp.logg.info(
                f"Time needed for {f.__name__}:"
                f" {time.strftime('%Hh:%Mm:%Ss', time.gmtime(end_time - start_time))}"
            )
        return ret

    return wrap


def get_file(key):
    """Resolve a package resource key to a filesystem path."""
    components = key.split("/")
    return str(
        importlib.resources.files(components[0]).joinpath("/".join(components[1:]))
    )


def with_timeout(timeout):
    """Return a decorator that limits a call to the given number of seconds."""

    def decorator(decorated):
        @functools.wraps(decorated)
        def inner(*args, **kwargs):
            with multiprocessing.pool.ThreadPool(1) as pool:
                async_result = pool.apply_async(decorated, args, kwargs)
                try:
                    return async_result.get(timeout)
                except multiprocessing.TimeoutError:
                    return None

        return inner

    return decorator


@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Update a tqdm progress bar as Joblib batches complete."""

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


def split_mut(mut):
    """Split an encoded mutation identifier into genomic components."""
    try:
        ref = mut.split(".")[-2]
        pos = mut.split(".")[-3]
        chrom = mut.split(".")[-4]
        gene = mut.split(".chr")[0].split("_")[1]
        ens = mut.split(".chr")[0].split("_")[0]
        alt = mut.split(".")[-1]
        return ens, gene, chrom, pos, ref, alt
    except IndexError:
        return None, None, None, None, None, None
