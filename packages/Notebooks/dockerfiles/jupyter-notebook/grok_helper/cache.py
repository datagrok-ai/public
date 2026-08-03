#!/usr/bin/env python

import os


class CacheEntity(object):
    def __init__(self, id='', version='0', start=0, end=0):
        self.id = id
        self.version = version
        self.start = start
        self.end = end

    def from_request(self, request):
        self.id = str(request.args.get('id', ''))
        self.version = request.args.get('version', '0', type=str)
        self.start = request.args.get('start', 0, type=int)
        self.end = request.args.get('end', 0, type=int)
        return self

    def from_json(self, map):
        self.id = '' if map['id'] is None else map['id']
        self.version = '0' if map['version'] is None else map['version']
        self.start = 0 if map['start'] is None else map['start']
        self.end = 0 if map['end'] is None else map['end']
        return self


class Cache(object):
    # Path to cache directory
    path = 'tmp'

    def __init__(self, path='tmp'):
        if not os.path.exists(path):
            os.mkdir(path)
        self.path = path

    def put(self, bytes, cache_entity) -> str:
        file_path = self.get_entity_path(cache_entity)
        file = open(file_path, 'wb+')
        file.write(bytes)
        file.close()
        return file_path

    def get(self, cache_entity):
        file = open(self.get_entity_path(cache_entity), 'rb')
        bytes = file.read()
        file.close()
        return bytes

    def get_column(self, cache_entity) -> list:
        file = open(self.get_entity_path(cache_entity), 'r')
        column = file.readlines()
        file.close()
        return column

    def exists(self, cache_entity):
        return os.path.isfile(self.get_entity_path(cache_entity))

    def remove(self, cache_entity):
        os.remove(self.get_entity_path(cache_entity))

    def clear(self, before_date):
        """
        Clears all entities before specified time

        :param before_date: Time, in milliseconds since epoch
        """
        before_date /= 1000.0
        for file in os.listdir(self.path):
            if os.path.getmtime(file) < before_date:
                os.remove(file)

    def get_entity_path(self, cache_entity, ext='.bin'):
        return os.path.join(self.path, cache_entity.id + '_' +
                            cache_entity.version + '_' +
                            str(cache_entity.start) + '-' + str(cache_entity.end) +
                            ext)
