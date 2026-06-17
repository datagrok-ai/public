#!/usr/bin/env python

import os
import subprocess


def call_process(params: list, working_dir: str) -> str:
    """
    Calls process.

    :param params: List of process parameters.
    :param working_dir: Working directory.
    :return: Log.
    """
    process = subprocess.Popen(params, cwd=working_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    log = os.fsdecode(out) + os.fsdecode(err)
    return log
