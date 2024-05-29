#!/usr/bin/env python
import logging
import os
import yaml
from multiprocessing import cpu_count

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

class Settings(object):
    """ GrokCompute server settings. """

    _configuration = None

    # General
    host = '0.0.0.0'
    port = 5000
    keys = None
    num_cores = 4
    cache_path = 'tmp'
    resources_path = 'resources'
    application_root = None

    # Modeling
    chemprop_path = None

    @staticmethod
    def init(configuration_path: str):
        conf = yaml.safe_load(open(configuration_path, 'r'))
        Settings._configuration = conf

        cache_path = conf['general'].get('cache_path')
        Settings.cache_path = os.path.join(os.getcwd(), 'tmp' if cache_path is None else cache_path)
        Settings.resources_path = os.path.join(os.getcwd(), conf['general'].get('resources_path'))

        host = conf['general'].get('host')
        Settings.host = '0.0.0.0' if host is None else host

        port = conf['general'].get('port')
        Settings.port = 5000 if port is None else port

        application_root = conf['general'].get('application_root')
        Settings.application_root = None if application_root is None else application_root

        keys = conf['general'].get('ssl_keys')
        if keys is not None:
            Settings.keys = tuple(keys.split(" "))

        cpu_limit = get_cpu_limit()
        num_cores = os.environ.get('GROK_COMPUTE_NUM_CORES')
        if num_cores is None:
            num_cores = conf['general'].get('num_cores')
        Settings.num_cores = cpu_limit if num_cores is None else int(num_cores)
        if Settings.num_cores > cpu_limit:
            logger.warning(
                f"Set CPU cores num to {cpu_limit} instead of {Settings.num_cores} because of {cpu_limit} CPU limit")
            Settings.num_cores = cpu_limit

        Settings.chemprop_path = os.path.join(os.getcwd(), conf['modeling'].get('chemprop_path'))