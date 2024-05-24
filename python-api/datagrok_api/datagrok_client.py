import requests
import pandas as pd
import json
from io import StringIO

class DatagrokClient:
    def __init__(self, api_key, base_url):
        self.api_key = api_key
        self.base_url = base_url
        self.headers = {
            "Authorization": self.api_key
        }

    def _request(self, method, endpoint, content_type="application/json", **kwargs):
        url = f"{self.base_url}{endpoint}"
        headers = self.headers.copy()
        headers["Content-Type"] = content_type
        response = requests.request(method, url, headers=headers, **kwargs)
        response.raise_for_status()
        return response

    def download_table(self, name):
        endpoint = f"/public/v1/tables/{name}"
        response = self._request('GET', endpoint, content_type="text/plain")
        return pd.read_csv(StringIO(response.text))

    def upload_table(self, name, dataframe):
        endpoint = f"/public/v1/tables/{name}"
        csv_data = dataframe.to_csv(index=False)
        response = self._request('POST', endpoint, content_type="text/csv", data=csv_data)
        return response.json()

    def download_file(self, connector, path):
        endpoint = f"/public/v1/files/{connector}/{path}"
        response = self._request('GET', endpoint, content_type="application/octet-stream")
        return pd.read_csv(StringIO(response.content.decode('utf-8')))

    def upload_file(self, connector, path, file_path):
        endpoint = f"/public/v1/files/{connector}/{path}"
        with open(file_path, 'rb') as file:
            response = self._request('POST', endpoint, content_type="application/octet-stream", data=file)
        return response.json()

    def share_dashboard(self, id, groups, access="View"):
        endpoint = f"/public/v1/dashboards/{id}/shares"
        params = {
            'groups': groups,
            'access': access
        }
        response = self._request('GET', endpoint, params=params)
        return response.json()

    def create_dashboard(self, name, table_ids, layout_filename=None):
        endpoint = f"/public/v1/dashboards/{name}/{table_ids}"
        if layout_filename:
            with open(layout_filename, 'r') as layout_file:
                layout = json.load(layout_file)
            response = self._request('POST', endpoint, json=layout, content_type="application/json")
        else:
            response = self._request('POST', endpoint, content_type="application/json")
        return response.json()

    def call_function(self, name, invocation_parameters):
        endpoint = f"/public/v1/{name}/call"
        response = self._request('POST', endpoint, json=invocation_parameters, content_type="application/json")
        return response.json()