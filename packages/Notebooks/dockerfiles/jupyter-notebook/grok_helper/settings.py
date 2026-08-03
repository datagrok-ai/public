#!/usr/bin/env python

import yaml


class Settings(object):
    """ GrokCompute server settings. """

    _configuration = None

    # General
    host = 'localhost'
    port = 5005
    keys = None
    num_cores = 4
    cache_path = 'tmp'
    environments_path = 'environments'
    application_root = None
    notebook_enabled = False

    @staticmethod
    def init(configuration_path: str):
        conf = yaml.safe_load(open(configuration_path, 'r'))
        Settings._configuration = conf

        cache_path = conf['general'].get('cache_path')
        Settings.cache_path = 'tmp' if cache_path is None else cache_path

        environments_path = conf['general'].get('environments_path')
        Settings.environments_path = 'environments' if environments_path is None else environments_path

        host = conf['general'].get('host')
        Settings.host = 'localhost' if host is None else host

        port = conf['general'].get('port')
        Settings.port = 5005 if port is None else port

        application_root = conf['general'].get('application_root')
        Settings.application_root = None if application_root is None else application_root

        notebook_enabled = conf['general'].get('notebook_enabled')
        Settings.notebook_enabled = False if notebook_enabled is None else notebook_enabled

        keys = conf['general'].get('ssl_keys')
        if keys is not None:
            Settings.keys = tuple(keys.split(" "))
