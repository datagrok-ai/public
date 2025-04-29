from typing import List
import requests
import pandas as pd
import json
from io import StringIO

from .group import Group, GroupMembershipRequest 

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
        name = name.replace(':', '.')
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
        name = name.replace(':', '.')
        endpoint = f"/public/v1/tables/{name}"
        csv_data = dataframe.to_csv(index=False).encode("utf-8")
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
        connector = connector.replace(':', '.')
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
        connector = connector.replace(':', '.')
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
        id = id.replace(':', '.')
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
        name = name.replace(':', '.')
        endpoint = f"/public/v1/{name}/call"
        response = self._request('POST', endpoint, json=invocation_parameters, content_type="application/json")
        return response.json()

    def list_groups(self, smart_filter: str = None, include_personal: bool = False, include_members: bool = False, include_memberships: bool = False) -> List["Group"]:
        '''
        Lists user groups from Datagrok.

        Parameters
        ----------
        smart_filter: str
            Optional. Smart search filter.
        include_personal: bool
            Optional. Include also personal groups. 
        include_members: bool
            Optional. Include group members. 
        include_memberships: bool
            Optional. Include group memberships.          

        Returns
        -------
        List[Group]
            A list of Group instances.
        '''

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
    
    def lookup_groups(self, query: str) -> List["Group"]:
        '''
        Looks up groups from Datagrok by name.

        Parameters
        ----------
        query: str
            Query string to match groups by name.

        Returns
        -------
        List[Group]
            A list of Group instances that match the query.
        '''
        params = {"query": query}
        endpoint = "/public/v1/groups/lookup"
        response = self._request('GET', endpoint, params=params)

        return [Group(self, **group_data) for group_data in response.json()]

    def get_group(self, group_id: str) -> "Group":
        '''
        Retrieves the details of a specific group by its ID. Includes group members and memberships.

        Parameters
        ----------
        group_id: str
            The ID of the group to retrieve details for.

        Returns
        -------
        Group
            A Group instance containing the details of the group.
        '''
        endpoint = f"/public/v1/groups/{group_id}"
        response = self._request('GET', endpoint)
        return Group(self, **response.json())
    
    def get_group_members(self, group_id: str, admin: bool = None) -> List["Group"]:
        '''
        Retrieves the members of a specific group by its ID.

        Parameters
        ----------
        group_id: str
            The ID of the group to retrieve members for.
        admin: bool, optional
            If set to True, only admin members will be included. If False, only non-admin members will be included.
            If not provided, only non-admin members are included.

        Returns
        -------
        List[Group]
            A list of Group instances representing the members of the group.
        '''

        endpoint = f"/public/v1/groups/{group_id}/members"
        params = {}
        if admin is not None:
            params["admin"] = admin
        response = self._request('GET', endpoint, params=params)
        return [Group(self, **group_data) for group_data in response.json()]

    def get_group_memberships(self, group_id: str, admin: bool = None) -> List["Group"]:
        '''
        Retrieves the memberships of a specific group by its ID.

        Parameters
        ----------
        group_id: str
            The ID of the group to retrieve members for.
        admin: bool, optional
            If set to True, only admin memberships will be included. If False, only non-admin memberships will be included.
            If not provided, only non-admin memberships are included.

        Returns
        -------
        List[Group]
            A list of Group instances representing the members of the group.
        '''
        # Prepare the URL and endpoint with the group ID
        endpoint = f"/public/v1/groups/{group_id}/memberships"
        
        # Prepare the query parameters
        params = {}
        if admin is not None:
            params["admin"] = admin

        response = self._request('GET', endpoint, params=params)
        return [Group(self, **group_data) for group_data in response.json()]

    def request_group_membership(self, group_id: str) -> GroupMembershipRequest:
        '''
        Submits a request to join a specific group by its ID.

        Parameters
        ----------
        group_id: str
            The ID of the group to request membership for.

        Returns
        -------
        None
            This method does not return any data but raises an exception if the request fails.
        '''

        endpoint = f"/public/v1/groups/{group_id}/request"
        response = self._request('POST', endpoint)
        return GroupMembershipRequest(self, **response.json())

    def get_group_membership_requests(self, group_id: str) -> List["GroupMembershipRequest"]:
        '''
        Retrieves a list of membership requests for a specific group.

        Parameters
        ----------
        group_id: str
            The ID of the group to retrieve membership requests for.

        Returns
        -------
        List[GroupMembershipRequest]
            A list of GroupMembershipRequest instances representing the requests for the group.
        '''

        endpoint = f"/public/v1/groups/{group_id}/requests"
        response = self._request('GET', endpoint)
        return [GroupMembershipRequest(self, **request_data) for request_data in response.json()]
    
    def approve_membership_request(self, request_id: str) -> None:
        """
        Approve a membership request for a group.

        Parameters
        ----------
        request_id: str
            The ID of the membership request to approve.

        Returns
        -------
        None
            This method performs an approval action and does not return anything.
        """
        endpoint = f"/public/v1/groups/requests/{request_id}/approve"
        self._request('POST', endpoint)
    
    def deny_membership_request(self, request_id: str) -> None:
        """
        Deny a membership request for a group.

        Parameters
        ----------
        request_id: str
            The ID of the membership request to deny.

        Returns
        -------
        None
            This method performs a denial action and does not return anything.
        """
        endpoint = f"/public/v1/groups/requests/{request_id}/deny"
        self._request('POST', endpoint)

    def get_current_user_group(self) -> Group:
        """
        Get the current user's group. Includes members and memberships

        Returns
        -------
        Group
            The current user's group instance.
        """
        endpoint = "/public/v1/user/group"
        response = self._request('GET', endpoint)
        return Group(self, **response.json())

    def save_group(self, group: Group, save_relations = False) -> Group:
        """
        Saves or updates group.

        Returns
        -------
        Group
            Created group
        """

        group.ensure_id() 
        endpoint = "/public/v1/groups"
        response = self._request('POST', endpoint, json=group.to_dict(), params={'saveRelations': str(save_relations).lower()})
        return Group(self, **response.json())
    
    def _request(self, method, endpoint, content_type="application/json", **kwargs):
        url = f"{self.base_url}{endpoint}"
        headers = self.headers.copy()
        headers["Content-Type"] = content_type
        response = requests.request(method, url, headers=headers, **kwargs)
        response.raise_for_status()
        content = response.content.decode()
        if 'ApiError' in content:
            raise Exception(content)
        return response
