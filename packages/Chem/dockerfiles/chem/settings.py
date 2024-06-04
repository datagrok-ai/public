import os
from multiprocessing import cpu_count
import logging


logger = logging.getLogger('grok_compute.settings')


def get_cpu_limit():
    try:
        with open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us") as fp:
            cfs_quota_us = int(fp.read())
        with open("/sys/fs/cgroup/cpu/cpu.cfs_period_us") as fp:
            cfs_period_us = int(fp.read())
        container_cpus = cfs_quota_us // cfs_period_us
    except FileNotFoundError as e:
        logger.warning('Can not get CPU cgroup files with restrictions')
        container_cpus = 0
    # For physical machine, the `cfs_quota_us` could be '-1'
    cpus = cpu_count() if container_cpus < 1 else container_cpus
    return cpus


class Settings:
    num_cores = 4
    port = 5000
    host = 'localhost'

    @staticmethod
    def init():
        host = os.environ['CHEM_HOST']
        Settings.host = 'localhost' if host is None else host
        port = os.environ['CHEM_PORT']
        Settings.port = 5000 if port is None else int(port)
        cpu_limit = get_cpu_limit()
        num_cores = os.environ.get('CHEM_NUM_CORES')
        Settings.num_cores = cpu_limit if num_cores is None else int(num_cores)
        print(Settings.num_cores)
        if Settings.num_cores > cpu_limit:
            logger.warning(
                f"Set CPU cores num to {cpu_limit} instead of {Settings.num_cores} because of {cpu_limit} CPU limit")
            Settings.num_cores = cpu_limit
