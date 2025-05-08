from typing import Optional, List, TYPE_CHECKING
from datetime import datetime

from datagrok_api.model import Model

if TYPE_CHECKING:
    from datagrok_client import DatagrokClient 

from typing import Optional, List
from datetime import datetime

class Group(Model):
    """Represents a group in the Datagrok system.
    
    A group can be a collection of other groups, with hierarchical relationships
    defined through parent-child relationships. Groups can have members and memberships,
    and can be personal.
    
    Attributes
    ----------
    name : Optional[str]
        The unique identifier name of the group
    friendly_name : Optional[str]
        A human-readable name for the group
    description : Optional[str]
        A description of the group's purpose
    created_on : Optional[datetime]
        When the group was created
    updated_on : Optional[datetime]
        When the group was last updated
    personal : bool
        Whether this is a personal group
    hidden : bool
        Whether this group is hidden from regular views
    parents : List[GroupRelation]
        The parent group relationships
    children : List[GroupRelation]
        The child group relationships
    """
    
    def __init__(self, client: Optional["DatagrokClient"] = None, **data):
        """Initialize a new Group instance.
        
        Parameters
        ----------
        client : Optional[DatagrokClient]
            The Datagrok client instance for making API calls
        **data : dict
            Dictionary containing group data including id, name, friendlyName, etc.
        """
        super().__init__(data.get("id"))
        self._client: "DatagrokClient" = client
        self.name: Optional[str] = data.get("name")
        self.friendly_name: Optional[str] = data.get("friendlyName")
        self.description: Optional[str] = data.get("description")

        created_str = data.get("createdOn")
        self.created_on: Optional[datetime] = (
            datetime.fromisoformat(created_str) if created_str else None
        )

        updated_str = data.get("updatedOn")
        self.updated_on: Optional[datetime] = (
            datetime.fromisoformat(updated_str) if updated_str else None
        )

        self.personal: bool = data.get("personal", False)
        self.hidden: bool = data.get("hidden", False)

        if data.get('parents') is None:
            data["parents"] = []

        if data.get('children') is None:
            data["children"] = []

        self.parents: List["GroupRelation"] = [
            GroupRelation(client, **rel) for rel in data.get("parents")
        ]

        self.children: List["GroupRelation"] = [
            GroupRelation(client, **rel) for rel in data.get("children")
        ]

    def _extract_groups(self, relations: List["GroupRelation"], admin: bool, direction: str) -> List["Group"]:
        """Extract groups from relations based on admin status and direction.
        
        Parameters
        ----------
        relations : List[GroupRelation]
            List of group relations to extract from
        admin : bool
            Whether to extract admin or non-admin groups
        direction : str
            Direction to extract groups from ('parent' or 'child')
            
        Returns
        -------
        List[Group]
            List of extracted groups
        """
        return [
            getattr(rel, direction)
            for rel in relations
            if rel.is_admin == admin and getattr(rel, direction) is not None
        ]

    @property
    def members(self) -> List["Group"]:
        """Get all non-admin members of this group.
        
        Returns
        -------
        List[Group]
            List of non-admin member groups
        """
        return self._extract_groups(self.children, admin=False, direction="child")

    @property
    def admin_members(self) -> List["Group"]:
        """Get all admin members of this group.
        
        Returns
        -------
        List[Group]
            List of admin member groups
        """
        return self._extract_groups(self.children, admin=True, direction="child")

    @property
    def memberships(self) -> List["Group"]:
        """Get all non-admin memberships of this group.
        
        Returns
        -------
        List[Group]
            List of non-admin membership groups
        """
        return self._extract_groups(self.parents, admin=False, direction="parent")

    @property
    def admin_memberships(self) -> List["Group"]:
        """Get all admin memberships of this group.
        
        Returns
        -------
        List[Group]
            List of admin membership groups
        """
        return self._extract_groups(self.parents, admin=True, direction="parent")
    
    def request_membership(self) -> "GroupMembershipRequest":
        """Request membership in this group.
        
        Returns
        -------
        GroupMembershipRequest
            The created membership request
        """
        return self._client.request_group_membership(self.id)
    
    def get_membership_requests(self) -> List["GroupMembershipRequest"]:
        """Get all pending membership requests for this group.
        
        Returns
        -------
        List[GroupMembershipRequest]
            List of pending membership requests
        """
        return self._client.get_group_membership_requests(self.id)
    
    def add_member(self, group: "Group", is_admin = False):
        """Add a member to this group.
        
        Parameters
        ----------
        group : Group
            The group to add as a member
        is_admin : bool, optional
            Whether the new member should have admin privileges, by default False
        """
        group.ensure_id()
        existing = next((c for c in self.children if c.child and c.child.id == group.id), None)

        if existing is not None:
            existing.is_admin = is_admin
        else:
            rel = GroupRelation(self._client, isAdmin=is_admin)
            rel.parent = Group(id=self.id)
            rel.child = group
            rel.ensure_id()
            self.children.append(rel)

        self._client.save_group(self, save_relations=True)        
    
    def to_dict(self):
        """Convert the group to a dictionary representation.
        
        Returns
        -------
        dict
            Dictionary containing all group data
        """
        return {
            "id": self.id,
            "name": self.name,
            "friendlyName": self.friendly_name,
            "description": self.description,
            "createdOn": self.created_on.isoformat() if self.created_on else None,
            "updatedOn": self.updated_on.isoformat() if self.updated_on else None,
            "personal": self.personal,
            "hidden": self.hidden,
            "parents": [rel.to_dict() for rel in self.parents],
            "children": [rel.to_dict() for rel in self.children],
        }

class GroupRelation(Model):
    """Represents a relationship between two groups.
    
    This class models the parent-child relationship between groups, including
    whether the relationship includes admin privileges.
    
    Attributes
    ----------
    is_admin : bool
        Whether this relationship includes admin privileges
    parent : Optional[Group]
        The parent group in the relationship
    child : Optional[Group]
        The child group in the relationship
    """
    
    def __init__(self, client: "DatagrokClient", **data):
        """Initialize a new GroupRelation instance.
        
        Parameters
        ----------
        client : DatagrokClient
            The Datagrok client instance for making API calls
        **data : dict
            Dictionary containing relation data including id, isAdmin, parent, and child
        """
        super().__init__(data.get("id"))
        self.is_admin: bool = data.get("isAdmin", False)

        self.parent: Optional[Group] = (
            Group(client, **data.get("parent")) if data.get("parent") else None
        )
        self.child: Optional[Group] = (
            Group(client, **data.get("child")) if data.get("child") else None
        )

    def to_dict(self):
        """Convert the group relation to a dictionary representation.
        
        Returns
        -------
        dict
            Dictionary containing all relation data
        """
        return {
            "id": self.id,
            "isAdmin": self.is_admin,
            "parent": self.parent.to_dict() if self.parent else None,
            "child": self.child.to_dict() if self.child else None,
        }    

class GroupMembershipRequest(Model):
    """Represents a request for group membership.
    
    This class models a request for a group to join another group, including
    the approval status and relevant timestamps.
    
    Attributes
    ----------
    approved : Optional[bool]
        Whether the request has been approved (None if pending)
    created_on : Optional[datetime]
        When the request was created
    updated_on : Optional[datetime]
        When the request was last updated
    resolution_date : Optional[datetime]
        When the request was resolved (approved/denied)
    from_group : Optional[Group]
        The group requesting membership
    to : Optional[Group]
        The group being requested to join
    """
    
    def __init__(self, client: "DatagrokClient", **data):
        """Initialize a new GroupMembershipRequest instance.
        
        Parameters
        ----------
        client : DatagrokClient
            The Datagrok client instance for making API calls
        **data : dict
            Dictionary containing request data including id, approved, timestamps, etc.
        """
        super().__init__(data.get("id"))
        self._client: "DatagrokClient" = client
        self.approved: Optional[bool] = data.get("approved")
        created_str = data.get("createdOn")
        self.created_on: Optional[datetime] = (
            datetime.fromisoformat(created_str) if created_str else None
        )

        updated_str = data.get("updatedOn")
        self.updated_on: Optional[datetime] = (
            datetime.fromisoformat(updated_str) if updated_str else None
        )

        resolution_str = data.get("resolutionDate")
        self.resolution_date: Optional[datetime] = (
            datetime.fromisoformat(resolution_str) if resolution_str else None
        )

        self.from_group: Optional[Group] = (
            Group(client, **data.get("from")) if data.get("from") else None
        )
        self.to: Optional[Group] = (
            Group(client, **data.get("to")) if data.get("to") else None
        )

    @property
    def is_viewed(self) -> bool:
        """Check if the request has been viewed (approved or denied).
        
        Returns
        -------
        bool
            True if the request has been viewed, False otherwise
        """
        return self.approved is not None
    
    def approve(self):
        """Approve this membership request."""
        return self._client.approve_membership_request(self.id)
    
    def deny(self):
        """Deny this membership request."""
        return self._client.deny_membership_request(self.id)
    
    def to_dict(self):
        """Convert the membership request to a dictionary representation.
        
        Returns
        -------
        dict
            Dictionary containing all request data
        """
        return {
            "id": self.id,
            "approved": self.approved,
            "createdOn": self.created_on.isoformat() if self.created_on else None,
            "updatedOn": self.updated_on.isoformat() if self.updated_on else None,
            "resolutionDate": self.resolution_date.isoformat() if self.resolution_date else None,
            "from": self.from_group.to_dict() if self.from_group else None,
            "to": self.to.to_dict() if self.to else None,
        }
