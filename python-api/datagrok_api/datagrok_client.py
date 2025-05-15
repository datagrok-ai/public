from typing import List, Optional, Dict, Any, Set, Union
import requests
import pandas as pd
import json
from io import StringIO
from pathlib import Path
import json
from requests_toolbelt.multipart import decoder
import os
import shutil

from .group import Group, GroupMembershipRequest 

class DatagrokClient:
    """A client for interacting with the Datagrok API.
    
    This class provides methods to interact with various Datagrok services including:
    - Table operations (download/upload)
    - File operations (download/upload/sync)
    - Dashboard management
    - Group management
    - Function calls
    
    The client handles authentication and provides a convenient interface for making
    API calls to the Datagrok server.
    
    Attributes
    ----------
    api_key : str
        The API key used for authentication. Should start with "Bearer".
    base_url : str
        The base URL of the Datagrok API (e.g., "https://public.datagrok.ai/api")
    headers : Dict[str, str]
        HTTP headers used for authentication in API requests
    """
    
    def __init__(self, api_key: str, base_url: str) -> None:
        """Initialize a new DatagrokClient instance.
        
        Parameters
        ----------
        api_key : str
            API key for authentication. Can be obtained from user profile page.
            Should start with "Bearer".
        base_url : str
            Base URL of the Datagrok API (e.g., "https://public.datagrok.ai/api")
        """
        self.api_key = api_key
        self.base_url = base_url
        self.headers = {
            "Authorization": self.api_key
        }

    def download_table(self, name: str) -> pd.DataFrame:
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
        response = self._request('GET', endpoint, content_type="text/plain")
        return pd.read_csv(StringIO(response.text))

    def upload_table(self, name: str, dataframe: pd.DataFrame) -> Dict[str, Any]:
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
        response = self._request('POST', endpoint, content_type="text/csv", data=csv_data)
        return response.json()

    def download_file(self, connector: str, path: str) -> Union[pd.DataFrame, bytes]:
        '''Download a file from Datagrok.

        Parameters
        ----------
        connector: str
            Identifier of a connector: ID or Grok name.
        path: str
            Path to file in a connector

        Returns
        -------
        Union[pd.DataFrame, bytes]
            If requested file is csv, returns corresponding pd.DataFrame.
            Otherwise returns bytes with file content.
        '''
        connector = connector.replace(':', '.')
        endpoint = f"/public/v1/files/{connector}/{path}"
        response = self._request('GET', endpoint, content_type="application/octet-stream")
        if path.endswith('.csv'):
            return pd.read_csv(StringIO(response.content.decode('utf-8')))
        return response.content

    def upload_file(self, connector: str, path: str, file_path: str) -> None:
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
        connector = connector.replace(':', '.')
        endpoint = f"/public/v1/files/{connector}/{path}"
        with open(file_path, 'rb') as file:
            self._request('POST', endpoint, content_type="application/octet-stream", data=file)

    def sync_dir(self, connector: str, remote_path: str, dir: str, ext: str = None, recursive: bool = True, sync_meta_file: str = ".sync_meta.json") -> None:
        """
        Syncs the content of a local directory with a remote directory on a server, updating local files and metadata based on the response.

        This method compares the current state of the local directory with the remote directory and updates files based on their 
        modified status or whether they have been deleted. It uses a provided connector to connect to a remote service and syncs the 
        data, storing metadata about file versions and deletions in a `.sync_meta.json` file. If a file has been deleted on the server, 
        it will also be deleted locally.

        Parameters
        ----------
        connector: str
            Identifier of a connector: ID or Grok name.

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
           

    def share_dashboard(self, id: str, groups: str, access: str = "View") -> None:
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
        id = id.replace(':', '.')
        endpoint = f"/public/v1/dashboards/{id}/shares"
        params = {
            'groups': groups,
            'access': access
        }
        self._request('GET', endpoint, params=params)

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
            response = self._request('POST', endpoint, json=layout, content_type="application/json")
        else:
            response = self._request('POST', endpoint, content_type="application/json")
        return response.json()

    def call_function(self, name: str, invocation_parameters: Dict[str, Any]) -> Any:
        """Call a Datagrok function with the specified parameters.
        
        Parameters
        ----------
        name : str
            Name of the function to call
        invocation_parameters : Dict[str, Any]
            Dictionary of parameters to pass to the function
            
        Returns
        -------
        Any
            The result of the function call
        """
        name = name.replace(':', '.')
        endpoint = f"/public/v1/{name}/call"
        response = self._request('POST', endpoint, json=invocation_parameters, content_type="application/json")
        return response.json()

    def list_groups(self, smart_filter: Optional[str] = None, include_personal: bool = False, 
                   include_members: bool = False, include_memberships: bool = False) -> List[Group]:
        """List groups from Datagrok with optional filtering and inclusion of related data.
        
        Parameters
        ----------
        smart_filter : Optional[str]
            Optional smart search filter to apply
        include_personal : bool, default=False
            Whether to include personal groups in the results
        include_members : bool, default=False
            Whether to include group members in the results
        include_memberships : bool, default=False
            Whether to include group memberships in the results
            
        Returns
        -------
        List[Group]
            List of Group instances matching the criteria
        """
        personal = str(include_personal).lower()
        smart_filter = smart_filter or f"personal={personal}"
        if smart_filter and not smart_filter.endswith(f"personal={personal}"):
            smart_filter += f" and personal={personal}"

        params = {"text": smart_filter}

        include_parts = []
        if include_members:
            include_parts.append("children.child")
        if include_memberships:
            include_parts.append("parents.parent")

        if include_parts:
            params["include"] = ",".join(include_parts)

        endpoint = "/public/v1/groups"
        response = self._request('GET', endpoint, params=params)
        
        return [Group(self, **group_data) for group_data in response.json()]
    
    def lookup_groups(self, query: str) -> List[Group]:
        """Search for groups by name.
        
        Parameters
        ----------
        query : str
            Search query to match against group names
            
        Returns
        -------
        List[Group]
            List of Group instances matching the search query
        """
        params = {"query": query}
        endpoint = "/public/v1/groups/lookup"
        response = self._request('GET', endpoint, params=params)

        return [Group(self, **group_data) for group_data in response.json()]

    def get_group(self, group_id: str) -> Group:
        """Get detailed information about a specific group.
        
        Parameters
        ----------
        group_id : str
            ID of the group to retrieve
            
        Returns
        -------
        Group
            Group instance with detailed information
        """
        endpoint = f"/public/v1/groups/{group_id}"
        response = self._request('GET', endpoint)
        return Group(self, **response.json())
    
    def get_group_members(self, group_id: str, admin: Optional[bool] = None) -> List[Group]:
        """Get members of a specific group.
        
        Parameters
        ----------
        group_id : str
            ID of the group to get members for
        admin : Optional[bool]
            If True, returns only admin members. If False, returns only non-admin members.
            If None, returns all members.
            
        Returns
        -------
        List[Group]
            List of Group instances representing the members
        """
        endpoint = f"/public/v1/groups/{group_id}/members"
        params = {}
        if admin is not None:
            params["admin"] = admin
        response = self._request('GET', endpoint, params=params)
        return [Group(self, **group_data) for group_data in response.json()]

    def get_group_memberships(self, group_id: str, admin: Optional[bool] = None) -> List[Group]:
        """Get memberships of a specific group.
        
        Parameters
        ----------
        group_id : str
            ID of the group to get memberships for
        admin : Optional[bool]
            If True, returns only admin memberships. If False, returns only non-admin memberships.
            If None, returns all memberships.
            
        Returns
        -------
        List[Group]
            List of Group instances representing the memberships
        """
        endpoint = f"/public/v1/groups/{group_id}/memberships"
        params = {}
        if admin is not None:
            params["admin"] = admin
        response = self._request('GET', endpoint, params=params)
        return [Group(self, **group_data) for group_data in response.json()]

    def request_group_membership(self, group_id: str) -> GroupMembershipRequest:
        """Request membership in a group.
        
        Parameters
        ----------
        group_id : str
            ID of the group to request membership for
            
        Returns
        -------
        GroupMembershipRequest
            The created membership request
        """
        endpoint = f"/public/v1/groups/{group_id}/request"
        response = self._request('POST', endpoint)
        return GroupMembershipRequest(self, **response.json())

    def get_group_membership_requests(self, group_id: str) -> List[GroupMembershipRequest]:
        """Get pending membership requests for a group.
        
        Parameters
        ----------
        group_id : str
            ID of the group to get membership requests for
            
        Returns
        -------
        List[GroupMembershipRequest]
            List of pending membership requests
        """
        endpoint = f"/public/v1/groups/{group_id}/requests"
        response = self._request('GET', endpoint)
        return [GroupMembershipRequest(self, **request_data) for request_data in response.json()]
    
    def approve_membership_request(self, request_id: str) -> None:
        """Approve a membership request.
        
        Parameters
        ----------
        request_id : str
            ID of the membership request to approve
        """
        endpoint = f"/public/v1/groups/requests/{request_id}/approve"
        self._request('POST', endpoint)
    
    def deny_membership_request(self, request_id: str) -> None:
        """Deny a membership request.
        
        Parameters
        ----------
        request_id : str
            ID of the membership request to deny
        """
        endpoint = f"/public/v1/groups/requests/{request_id}/deny"
        self._request('POST', endpoint)

    def get_current_user_group(self) -> Group:
        """Get the current user's group.
        
        Returns
        -------
        Group
            Group instance representing the current user's group
        """
        endpoint = "/public/v1/user/group"
        response = self._request('GET', endpoint)
        return Group(self, **response.json())

    def save_group(self, group: Group, save_relations: bool = False) -> Group:
        """Save or update a group.
        
        Parameters
        ----------
        group : Group
            The group to save or update
        save_relations : bool, default=False
            Whether to save group relations (members and memberships)
            
        Returns
        -------
        Group
            The saved group instance
        """
        group.ensure_id() 
        endpoint = "/public/v1/groups"
        response = self._request('POST', endpoint, json=group.to_dict(), params={'saveRelations': str(save_relations).lower()})
        return Group(self, **response.json())
    
    def _request(self, method: str, endpoint: str, content_type: str = "application/json", **kwargs: Any) -> requests.Response:
        """Make an HTTP request to the Datagrok API.
        
        Parameters
        ----------
        method : str
            HTTP method to use (GET, POST, etc.)
        endpoint : str
            API endpoint to call
        content_type : str, default="application/json"
            Content type for the request
        **kwargs : Any
            Additional arguments to pass to requests.request
            
        Returns
        -------
        requests.Response
            The response from the API
            
        Raises
        ------
        Exception
            If the API returns an error or if the response is not successful
        """
        url = f"{self.base_url}{endpoint}"
        headers = self.headers.copy()
        headers["Content-Type"] = content_type
        response = requests.request(method, url, headers=headers, **kwargs)
        response.raise_for_status()
        if 'api-error' in response.headers:
            content = response.content.decode()
            if 'ApiError' in content:
                raise Exception(content)
            else:
                raise Exception("Something went wrong during request") 
        return response
    
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

