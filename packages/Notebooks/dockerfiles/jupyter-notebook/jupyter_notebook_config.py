# Configuration file for jupyter-notebook.

import os
import io
import json
import urllib3
import requests

c.NotebookApp.base_url = '/notebook'
c.NotebookApp.allow_credentials = True
c.NotebookApp.allow_origin = '*'
c.NotebookApp.allow_root = True
c.NotebookApp.disable_check_xsrf = True
c.NotebookApp.ip = '0.0.0.0'
c.NotebookApp.notebook_dir = '/home/grok/notebooks'
c.NotebookApp.port = 8889
c.NotebookApp.token = '' #4b8ed936cf61a1c5e37a8b3a845599941c272de6e29330a0
c.NotebookApp.password = ''
c.NotebookApp.tornado_settings = {'headers': {'Content-Security-Policy': 'frame-ancestors *'}}
c.MappingKernelManager.cull_idle_timeout = 600
c.NotebookApp.open_browser = False

def script_post_save(model, os_path, contents_manager, **kwargs):
    if model['type'] != 'notebook':
        return

    nb_file = io.open(os_path, 'r')
    nb_str = nb_file.read()
    nb_file.close()
    nb = json.loads(nb_str)
    if 'datagrok' in nb['metadata']:
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        params = nb['metadata']['datagrok']
        host = params['host']
        if 'localhost' in host and os.environ.get('CVM_ENVIRONMENT') == 'docker':
            host = host.replace('localhost', 'host.docker.internal')
        token = params['session_token']
        requests.post(host + "/notebooks/file/" + params['id'],
                      data=nb_str.encode('utf-8'),
                      headers={
                          "Authorization": token
                      }, verify=False)


c.FileContentsManager.post_save_hook = script_post_save
