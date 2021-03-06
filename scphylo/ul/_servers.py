import csv
import os

import scphylo as scp


def write_cmds_get_main(
    df, jobname, time, mem, module, nthread, email, tmpdir, afterok=None
):
    cmdfile = f"{tmpdir}/cmd/{jobname}.sh"
    logdir = f"{tmpdir}/log/{jobname}"
    scp.ul.mkdir(f"{tmpdir}/cmd")
    scp.ul.mkdir(logdir)
    if "biowulf" in os.uname()[1]:
        if df.shape[0] > 1:
            cmd = [
                "swarm",
                f"--file {cmdfile}",
                f"--time {time}",
                f"--gb-per-process {mem}",
                f"--logdir {logdir}",
                "--bundle 1",
                "--partition ccr,norm",
                "--processes-per-subjob 1",
                f"--threads-per-process {nthread}",
                "--noht",
            ]
            if module is not None:
                cmd += [
                    f"--module {module}",
                ]
            tmp = [f'--sbatch "--mail-type=END --mail-user={email}']
            if afterok is not None:
                tmp += [f"--dependency=afterok:{afterok}"]
            if jobname == "sra":
                tmp += ["--gres=lscratch:20"]
            tmp += [f'--job-name={jobname}"']
            cmd += tmp
        else:
            cmd = [
                "sbatch",
                "--ntasks=1",
                f"--cpus-per-task={nthread}",
                "--ntasks-per-core=1",
                f"--time={time}",
                f"--mem={mem}g",
                f"--job-name={jobname}",
                "--partition=ccr,norm",
                f"--chdir={tmpdir}",
                f"--output={logdir}/log.o",
                f"--error={logdir}/log.e",
                "--mail-type=END",
                f"--mail-user={email}",
            ]
            if afterok is not None:
                cmd += [f"--dependency=afterok:{afterok}", f"{cmdfile}"]
            else:
                cmd += [f"{cmdfile}"]
    else:
        cmd = [
            f"cat {cmdfile} | parallel -j0 --bar {{}}",
        ]
    cmdmain = " ".join(cmd)
    with open(cmdfile, "w") as out:
        out.write("#!/bin/bash\n\n#" + cmdmain + "\n\n")
        df["cmd"].to_csv(
            out, index=False, header=False, quoting=csv.QUOTE_NONE, sep="\t"
        )
    return cmdmain


def cmd(arr, islast=False):
    if islast:
        return " ".join(arr)
    else:
        return " ".join(arr) + " && "
