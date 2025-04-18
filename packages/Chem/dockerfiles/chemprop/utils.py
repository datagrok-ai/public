import multiprocessing
import subprocess
import logging
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)

class Types(object):
    """ Commonly used type identifiers. """
    INT = 'int'
    FLOAT = 'double'
    BOOL = 'bool'
    STRING = 'string'
    DATE_TIME = 'datetime'
    LIST = 'list'

def call_process(params: list):
    """
    Calls a process with the given parameters.

    :param params: List of process parameters.
    :return: Log of the process's stdout if successful.
    :raises RuntimeError: If the process fails (return code != 0).
    """
    logger.debug(f'Params: {params}')
    process = subprocess.run(params, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    output = process.stdout.decode()

    if process.returncode != 0:
        logger.error(f'Process failed with output: {output}')
        raise RuntimeError(f"Process failed with return code {process.returncode}: {output}")

    logger.debug(f'Process succeeded with output: {output}')
    return output
