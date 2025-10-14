import json
import os
import shutil
from io import StringIO
from pathlib import Path
from typing import List, Optional, Set, Union

import pandas as pd
from requests_toolbelt.multipart import decoder

from datagrok_api.http_client import HttpClient
from datagrok_api.models.data_connection import DataConnection, DataSourceType


class FilesClient:

    def __init__(self, client: HttpClient):
        self.client = client    

    def download(self, connector: Union[DataConnection, str], path: str) -> Union[pd.DataFrame, bytes]:
        '''Download a file from Datagrok.

        Parameters
        ----------
        connector: str or DataConnection
            DataConnection or ID or Grok name.
        path: str
            Path to file in a connector

        Returns
        -------
        Union[pd.DataFrame, bytes]
            If requested file is csv, returns corresponding pd.DataFrame.
            Otherwise returns bytes with file content.
        '''
        if isinstance(connector, DataConnection):
            self._validate_connection(connector)
            connector = connector.grok_name
        connector = connector.replace(':', '.')
        endpoint = f"/public/v1/files/{connector}/{path}"
        response = self.client.get(endpoint, headers={'Content-Type': 'application/octet-stream'})
        if path.endswith('.csv'):
            return pd.read_csv(StringIO(response.content.decode('utf-8')))
        return response.content
    
    def upload(self, connector: Union[DataConnection, str], path: str, file_path: str) -> None:
        '''
        Uploads a file to Datagrok.

        Parameters
        ----------
        connector: str or DataConnection
            DataConnection or ID or Grok name.
        path: str
            Path to file in a connector
        file_path: str
            Path to local file that needs to be uploaded
        '''
        if isinstance(connector, DataConnection):
            self._validate_connection(connector)
            connector = connector.grok_name
        connector = connector.replace(':', '.')
        endpoint = f"/public/v1/files/{connector}/{path}"
        with open(file_path, 'rb') as file:
            self.client.post(endpoint, headers={'Content-Type': 'application/octet-stream'}, data=file)

    def sync_dir(self, connector: Union[DataConnection, str], remote_path: str, dir: str, ext: str = None, recursive: bool = True, sync_meta_file: str = ".sync_meta.json") -> None:
        """
        Syncs the content of a local directory with a remote directory on a server, updating local files and metadata based on the response.

        This method compares the current state of the local directory with the remote directory and updates files based on their 
        modified status or whether they have been deleted. It uses a provided connector to connect to a remote service and syncs the 
        data, storing metadata about file versions and deletions in a `.sync_meta.json` file. If a file has been deleted on the server, 
        it will also be deleted locally.

        Parameters
        ----------
        connector: str or DataConnection
            DataConnection or ID or Grok name.

        remote_path: str
            The path on the server connection to sync with. Should be directory

        dir: str
            The local directory to sync. If the directory doesn't exist, it will be created.

        ext: str, optional
            A file extension filter to limit the types of files to sync. If not specified, all files are synced.

        recursive: bool
            Whether to sync the directory recursively, including subdirectories. If False, only files in the root of the directory are synced.

        sync_meta_file: str, optional, default=".sync_meta.json"
            The name of the metadata file that stores the eTags and sync information for the files. This file is used to track the current 
            state of the local directory and should be present for syncing to work correctly.
        """
        if isinstance(connector, DataConnection):
            self._validate_connection(connector)
            connector.ensure_id()
            connector = connector.id
            
        os.makedirs(dir, exist_ok=True)
        sync_meta_file_path = Path(dir, sync_meta_file)
        sync_meta_dict = {}
        if sync_meta_file_path.exists():
            with open(sync_meta_file_path, "r", encoding="utf-8") as f:
                content = f.read()
                sync_meta_dict = json.loads(content)

        existed = self._list_files_unix_style(dir, {sync_meta_file})
        for key in list(sync_meta_dict.keys()):
            if key not in existed:
                del sync_meta_dict[key]

        if remote_path.startswith('/'):
            remote_path = remote_path[1:]
            
        connector = connector.replace(':', '.')
        endpoint = f"/public/v1/files/sync/{connector}/{remote_path}"
        params = {'recursive': str(recursive).lower()}
        if ext:
           params['ext'] = ext 

        response = self._request('POST', endpoint, content_type="application/json", params=params, json=sync_meta_dict)

        def update_sync_meta(content: dict):
            with open (sync_meta_file_path, "w") as f:
                f.write(json.dumps(content))   


        # remote_path is empty folder    
        if response.status_code == 204:
            self._clear_directory(dir)
            update_sync_meta({})  
            return
        
        if not response.headers.get('Content-Type', '').startswith('multipart/'):
            raise Exception(f"Incorrect response type. Expected multipart but received {response.headers.get('Content-Type', '')}")
        
        multipart_data = decoder.MultipartDecoder.from_response(response)

        deleted_files = set()

        for i, part in enumerate(multipart_data.parts):
            content_disposition = part.headers.get(b'Content-Disposition', b'').decode()
            filename = f"part_{i}"
            if "filename=" in content_disposition:
                filename = content_disposition.split("filename=")[1].strip('"')

            status_bytes = part.headers.get(b'X-File-Status')
            status = status_bytes.decode('utf-8') if status_bytes else None

            if status == "deleted":
                deleted_files.add(filename)
                continue
            elif status == "not-modified":
                continue
       
            full_filepath = os.path.join(dir, filename.replace('/', os.sep))   
            os.makedirs(os.path.dirname(full_filepath), exist_ok=True)

            with open(full_filepath, "wb") as f:
                f.write(part.content)

            e_tag = part.headers.get(b'X-File-E-Tag')
            if e_tag:
                sync_meta_dict[filename] = e_tag.decode('utf-8')


        for filename in deleted_files:
            del sync_meta_dict[filename]
            full_filepath = Path(dir, filename.replace('/', os.sep))
            if full_filepath.exists():
                full_filepath.unlink()

        update_sync_meta(sync_meta_dict)

    def _list_files_unix_style(self, root_dir: str, excludes: Optional[Set] = None) -> List[str]:
        """List files in a directory using Unix-style paths.
        
        Parameters
        ----------
        root_dir : str
            Root directory to list files from
        excludes : set, optional
            Set of file paths to exclude from the listing
            
        Returns
        -------
        List[str]
            List of file paths relative to the root directory
        """
        root = Path(root_dir)

        if excludes is None:
            excludes = set()

        files = []
        for path in root.rglob('*'):
            if path.is_file():
                relative_path = path.relative_to(root).as_posix()
                if relative_path not in excludes:
                    files.append(relative_path)

        return files

    def _clear_directory(self, path: str):
        """Remove all files and subdirectories from a directory.
        
        Parameters
        ----------
        path : str
            Path to the directory to clear
        """
        dir_path = Path(path)
        for item in dir_path.iterdir():
            if item.is_file() or item.is_symlink():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)

    def _validate_connection(self, conn: DataConnection):
        if not DataSourceType.is_file_data_source(conn.data_source):
            raise Exception('Files operations are only available for File type data sources')             
