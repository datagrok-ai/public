import requests
import pandas as pd
import json
from io import StringIO

class DatagrokClient:
    def __init__(self, api_key, base_url):
        '''
        Create a datagrok api client.

        Parameters
        ----------
        api_key: str
            Api key for a client. Can be obtained from user profile page.
            Starts with "Bearer".
        base_url: str
            Url of datagrok API. For example, "https://public.datagrok.ai/api"
        
        Returns
        -------
        DatagrokClient
            New instance of Datagrok Client
        '''
        self.api_key = api_key
        self.base_url = base_url
        self.headers = {
            "Authorization": self.api_key
        }

    def download_table(self, name):
        '''
        Downloads a table from Datagrok.

        Parameters
        ----------
        name: str
            Identifier of a table. Can be accessed from context panel.
            Several options supppored: UUID, Grok names (e.g. Namespace.Project.Table)
        
        Returns
        -------
        pandas.DataFrame
            Dataframe with table data
        '''
        endpoint = f"/public/v1/tables/{name}"
        response = self._request('GET', endpoint, content_type="text/plain")
        return pd.read_csv(StringIO(response.text))

    def upload_table(self, name, dataframe):
        '''
        Uploads a table to Datagrok.

        Parameters
        ----------
        name: str
            Name for new table.
            If in grok format, can also specify project and namespace.
        dataframe: pandas.Dataframe
            table data to save

        Returns
        -------
        object
            set of identifiers for the uploaded table
        '''
        endpoint = f"/public/v1/tables/{name}"
        csv_data = dataframe.to_csv(index=False)
        response = self._request('POST', endpoint, content_type="text/csv", data=csv_data)
        return response.json()

    def download_file(self, connector, path):
        '''
        Download a file from Datagrok.

        Parameters
        ----------
        connector: str
            Identifier of a connector: ID or Grok name.
        path: str
            Path to file in a connector

        Returns
        -------
        object
            If requested file is csv, returns corresponding pd.Datafrme.
            Otherwise returns list of bytes with file content.
        '''
        endpoint = f"/public/v1/files/{connector}/{path}"
        response = self._request('GET', endpoint, content_type="application/octet-stream")
        if path.endswith('.csv'):
            return pd.read_csv(StringIO(response.content.decode('utf-8')))
        return response.content

    def upload_file(self, connector, path, file_path):
        '''
        Uploads a file to Datagrok.

        Parameters
        ----------
        connector: str
            Identifier of a connector: ID or Grok name.
        path: str
            Path to file in a connector
        file_path: str
            Path to local file that needs to be uploaded
        '''
        endpoint = f"/public/v1/files/{connector}/{path}"
        with open(file_path, 'rb') as file:
            self._request('POST', endpoint, content_type="application/octet-stream", data=file)

    def share_dashboard(self, id, groups, access="View"):
        '''
        Shares a dashboard.

        Parameters
        ----------
        id: str
            Identifier of a project: ID or Grok name.
        groups: str
            Comma-separated list of group/user names
        access: str
            Either 'View' or 'Edit'
        '''
        endpoint = f"/public/v1/dashboards/{id}/shares"
        params = {
            'groups': groups,
            'access': access
        }
        self._request('GET', endpoint, params=params)

    def create_dashboard(self, name, table_ids, layout_filename=None):
        '''
        Shares a dashboard.

        Parameters
        ----------
        name: str
            Name for new dashboard.
        table_ids: str
            Comma-separated list of table ids
        layout_filename: str
            Optional. Filename for a project layout.

        Returns
        -------
        object
            Identifiers of a project.
        '''
        endpoint = f"/public/v1/dashboards/{name}/{table_ids}"
        if layout_filename:
            with open(layout_filename, 'r') as layout_file:
                layout = json.load(layout_file)
            response = self._request('POST', endpoint, json=layout, content_type="application/json")
        else:
            response = self._request('POST', endpoint, content_type="application/json")
        return response.json()

    def call_function(self, name, invocation_parameters):
        '''
        Performs a call of datagrok function.

        Parameters
        ----------
        name: str
            Function name.
        invocation_parameters: dict
            Dict with parameter values for function call.

        Returns
        -------
        object
            Result of invocation: either a single value or a list or outputs.
        '''
        endpoint = f"/public/v1/{name}/call"
        response = self._request('POST', endpoint, json=invocation_parameters, content_type="application/json")
        return response.json()

    def _request(self, method, endpoint, content_type="application/json", **kwargs):
        url = f"{self.base_url}{endpoint}"
        headers = self.headers.copy()
        headers["Content-Type"] = content_type
        response = requests.request(method, url, headers=headers, **kwargs)
        response.raise_for_status()
        return response
