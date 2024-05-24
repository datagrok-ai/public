#!/usr/bin/env python

import io
import urllib3
import pandas as pd
import requests

grok_host = ''
grok_write_token = ''
grok_user_id = ''

def grok_init(host, write_token='', user_id=''):
    """
    Init Datagrok API module.

    :param host: Datagrok server URL.
    :param write_token: Table write token.
    :param user_id: User ID.
    """
    global grok_host, grok_write_token, grok_user_id
    grok_host = host
    grok_write_token = write_token
    grok_user_id = user_id


def grok_read(token):
    """
    Reads table from Datagrok.

    :param host: Datagrok server URL.
    :param token: Table token.
    :return: Pandas DataFrame.
    """
    global grok_host
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    session = requests.Session()
    session.verify = False
    response = session.get(f'{grok_host}/tables/csv/{token}.csv')
    return pd.read_csv(io.StringIO(response.content.decode('utf-8')))


def grok(table):
    """
    Writes table to Datagrok server.

    :param table: Pandas DataFrame.
    """
    global grok_host, grok_write_token, grok_user_id
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    csv = table.to_csv(index=False)
    url = f'{grok_host}/tables/csv/{grok_write_token}'
    if grok_user_id is not None:
        url = f'{url}?userId={grok_user_id}'
    requests.post(url, data=csv, verify=False)
