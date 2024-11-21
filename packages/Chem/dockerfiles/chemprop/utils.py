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
    process = subprocess.run(params, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if process.returncode != 0:
        error_log = process.stderr.decode()
        logger.error(f'Process failed with error log: {error_log}')
        raise RuntimeError(f"Process failed with return code {process.returncode}: {error_log}")
    
    log = process.stdout.decode()
    logger.debug(f'Process succeeded with log: {log}')
    return log
