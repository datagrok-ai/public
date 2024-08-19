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
    Calls process.

    :param params: List of process parameters.
    :param working_dir: Working directory.
    :return: Log.
    """
    process = subprocess.run(params, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    log = process.stderr
    logger.debug(f'Process log: {log}')
    return log