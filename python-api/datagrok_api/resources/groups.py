from typing import List, Optional, Union

from datagrok_api.http_client import HttpClient
from datagrok_api.models.group import Group
from datagrok_api.models.user import User


class GroupsClient:
    def __init__(self, client: HttpClient):
        self.client = client

    def find(self, query: str) -> List[Group]:
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
        response = self.client.get(endpoint, params=params)

        return [Group.from_dict(group_data) for group_data in response.json()]
    
    def get(self, id: str) -> Group:
        """Get detailed information about a specific group.
        
        Parameters
        ----------
        id : str
            ID or name of the group to retrieve
            
        Returns
        -------
        Group
            Group instance with detailed information
        """
        id = id.replace(':', '.')
        endpoint = f"/public/v1/groups/{id}"
        response = self.client.get(endpoint)
        return Group.from_dict(response.json())
    
    def save(self, group: Group, save_relations: bool = False) -> Group:
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
        response = self.client.post(endpoint, json=group.to_dict(), params={'saveRelations': str(save_relations).lower()})
        return Group.from_dict(response.json())

    def delete(self, group: Union[Group, str]):
        """Deletes a group.
        
        Parameters
        ----------
        group : Group or str
            The group to delete or (group name or group id)
        """
        if isinstance(group, Group):
            group.ensure_id()
            group = group.id
        group = group.replace(':', '.')    
        endpoint = f"/public/v1/groups/{group}"
        self.client.delete(endpoint)

    def list(self, smart_filter: Optional[str] = None, include_personal: bool = False, 
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
        response = self.client.get(endpoint, params=params)
        
        return [Group.from_dict(group_data) for group_data in response.json()]
    
    def add_member(self, parent: Union[Group, str], 
                   child: Union[Group, User, str], is_admin: bool = False) -> Group:
        """
        Adds a child group or user as a member of a parent group.

        Parameters
        ----------
        parent : Group or str
            The parent group or its name.
        child : Group, User, or str
            The group, user, or their name to be added to the parent group.
        is_admin : bool, optional
            Whether the child should be granted admin privileges (default is False).

        Returns
        -------
        Group
            The updated parent group with the new member added.
        """

        def resolve_group(name_or_group: Union[str, Group]) -> Group:
            if isinstance(name_or_group, Group):
                return name_or_group
            matches = self.find(name_or_group)
            if len(matches) != 1:
                raise Exception(f"Can't find group '{name_or_group}' or name is ambiguous")
            return matches[0]

        if isinstance(parent, str):
            parent = resolve_group(parent)

        if isinstance(child, User):
            matches = self.find(child.name)
            child_group = next((g for g in matches if g.personal), None)
            if not child_group:
                raise Exception(f"Can't find personal group for user '{child.name}'")
            child = child_group

        elif isinstance(child, str):
            child = resolve_group(child)

        elif not isinstance(child, Group):
            raise TypeError(f"Unsupported type for child: {type(child)}")

        parent.add_member(child, is_admin=is_admin)
        return self.save(parent, save_relations=True)

    def get_members(self, group: Group, admin: Optional[bool] = None) -> List[Group]:    
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
        endpoint = f"/public/v1/groups/{group.id}/members"
        params = {}
        if admin is not None:
            params["admin"] = admin
        response = self.client.get(endpoint, params=params)
        return [Group.from_dict(group_data) for group_data in response.json()]
    
    def get_memberships(self, group: Group, admin: Optional[bool] = None) -> List[Group]:
        """Get memberships of a specific group.
        
        Parameters
        ----------
        group : Group
            The group to get memberships for
        admin : Optional[bool]
            If True, returns only admin memberships. If False, returns only non-admin memberships.
            If None, returns all memberships.
            
        Returns
        -------
        List[Group]
            List of Group instances representing the memberships
        """
        endpoint = f"/public/v1/groups/{group.id}/memberships"
        params = {}
        if admin is not None:
            params["admin"] = admin
        response = self.client.get(endpoint, params=params)
        return [Group.from_dict(group_data) for group_data in response.json()]

    def current(self) -> Group:
        """Returns the group associated with the current authenticated user.
        
        Returns
        -------
        Group
            Group instance representing the current user's group
        """
        endpoint = "/public/v1/groups/current"
        response = self.client.get(endpoint)
        return Group.from_dict(response.json())
