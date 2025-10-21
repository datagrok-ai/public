import json
from io import StringIO
from typing import Any, Dict, Optional

import pandas as pd

from datagrok_api.http_client import HttpClient


class TablesClient:

    def __init__(self, client: HttpClient):
        self.client = client

    def download(self, name: str) -> pd.DataFrame:
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
        name = name.replace(':', '.')
        endpoint = f"/public/v1/tables/{name}"
        response = self.client.get(endpoint, headers={'Content-Type': 'text/plain'})
        return pd.read_csv(StringIO(response.text))
    
    def upload(self, name: str, dataframe: pd.DataFrame) -> Dict[str, Any]:
        '''
        Uploads a table to Datagrok.

        Parameters
        ----------
        name: str
            Name for new table. If in grok format, can also specify project and namespace.
        dataframe: pandas.Dataframe
            table data to save

        Returns
        -------
        Dict[str, Any]
            set of identifiers for the uploaded table
        '''
        name = name.replace(':', '.')
        endpoint = f"/public/v1/tables/{name}"
        csv_data = dataframe.to_csv(index=False).encode("utf-8")
        response = self.client.post(endpoint, headers={'Content-Type': 'text/csv'}, data=csv_data)
        return response.json()
    
    def create_dashboard(self, name: str, table_ids: str, layout_filename: Optional[str] = None) -> Dict[str, Any]:
        """Create a new dashboard in Datagrok.
        
        Parameters
        ----------
        name : str
            Name for the new dashboard
        table_ids : str
            Comma-separated list of table IDs to include in the dashboard
        layout_filename : Optional[str]
            Optional path to a JSON file containing the dashboard layout
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing the dashboard identifiers and metadata
        """
        name = name.replace(':', '.')
        table_ids = table_ids.replace(':', '.')
        endpoint = f"/public/v1/dashboards/{name}/{table_ids}"
        if layout_filename:
            with open(layout_filename, 'r') as layout_file:
                layout = json.load(layout_file)
            response = self.client.post(endpoint, json=layout)
        else:
            response = self.client.post(endpoint)
        return response.json()
