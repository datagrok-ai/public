from typing import Any, Dict, Optional, List
from datetime import datetime

from datagrok_api.models.model import Model, NamedModel

from typing import Optional, List
from datetime import datetime


class Group(NamedModel):
    """
    Represents a group in the Datagrok system.

    A group can be a collection of other groups, with hierarchical relationships
    defined through parent-child links. Groups can be personal, hidden, and have
    admin or non-admin members.

    Attributes
    ----------
    description : Optional[str]
        A description of the group's purpose.
    personal : bool, optional
        Whether this is a personal group of the user (default is False).
    hidden : bool, optional
        Whether this group is hidden from the UI (default is False).
        name : str
        Internal name of the model, used as an identifier.
    id : Optional[str]
        Unique identifier for the model instance. Can be None for new instances.
    friendly_name : Optional[str]
        Human-readable name for display purposes. Defaults to None.
    created_on : Optional[datetime]
        Timestamp indicating when the model was created. Defaults to None.
    updated_on : Optional[datetime]
        Timestamp indicating when the model was last updated. Defaults to None.
    namespace : Optional[str]
        Optional namespace prefix to distinguish models with the same name. Defaults to None.
    """

    _all_users: 'Group' = None
    _developers: 'Group' = None
    _admin: 'Group' = None
    _system: 'Group' = None
    _administrators: 'Group' = None
    _test: 'Group' = None

    def __init__(
        self,
        description: Optional[str] = None,
        hidden: bool = False,
        personal: bool = False,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.description = description
        self.hidden = hidden
        self.personal = personal
        self.parents: List["GroupRelation"] = []
        self.children: List["GroupRelation"] = []


    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Group":
        """
        Internal method: create a Group from a dictionary (e.g., API response).

        Parameters
        ----------
        client : Optional[DatagrokClient]
            The client used to fetch related entities.
        data : dict
            Group data in dictionary form.

        Returns
        -------
        Group
            The deserialized Group object.
        """
        group = cls(
            id=data.get("id"),
            name=data.get("name") or "",
            friendly_name=data.get("friendlyName"),
            description=data.get("description"),
            hidden=data.get("hidden", False),
            personal=data.get("personal", False),
            namespace=data.get("namespace")
        )

        created_str = data.get("createdOn")
        group.created_on = datetime.fromisoformat(created_str) if created_str else None

        updated_str = data.get("updatedOn")
        group.updated_on = datetime.fromisoformat(updated_str) if updated_str else None

        group.parents = [
            GroupRelation.from_dict(rel) for rel in data.get("parents") or []
        ]
        group.children = [
            GroupRelation.from_dict(rel) for rel in data.get("children") or []
        ]

        return group

    @classmethod
    def ALL_USERS(cls):
        if cls._all_users is None:
           cls._all_users = Group(id="a4b45840-9a50-11e6-9cc9-8546b8bf62e6", name="All users")
        return Group._all_users
    
    @classmethod
    def DEVELOPERS(cls):
        if cls._developers is None:
            cls._developers = Group(id="ba9cd191-9a50-11e6-9cc9-910bf827f0ab", name="Developers")
        return Group._developers
    
    @classmethod
    def TEST(cls):
        if cls._test is None:
            cls._test = Group(id="ca1e672e-e3be-40e0-b79b-8546b8bf62e6", name="Test")
        return Group._test
    
    @classmethod
    def ADMIN(cls):
        if cls._admin is None:
            cls._admin = Group(id="a4b45840-9a50-11e6-c537-6bf8e9ab02ee", name="Admin")
        return Group._admin
    
    @classmethod
    def SYSTEM(cls):
        if cls._system is None:
            cls._system = Group(id="a4b45840-ac9c-4d39-8b4b-4db3e576b3c3", name="System")
        return Group._system
    
    @classmethod
    def ADMINISTRATORS(cls):
        if cls._administrators is None:
            cls._administrators = Group(id="1ab8b38d-9c4e-4b1e-81c3-ae2bde3e12c5", name="Administrators")
        return Group._administrators

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
    
    def add_member(self, group: "Group", is_admin = False):
        """Adds a child group as a member of this group (in-memory only).
        Group should be save using client to actually add a member.
        
        Parameters
        ----------
        group : Group
            The group to add as a member
        is_admin : bool, optional
            Whether the new member should have admin privileges, by default False
        """
        self.ensure_id()
        group.ensure_id()
        existing = next((c for c in self.children if c.child and c.child.id == group.id), None)

        if existing is not None:
            existing.is_admin = is_admin
        else:
            rel = GroupRelation(parent=Group(id=self.id, name=self.name), child=group, is_admin=is_admin)
            rel.ensure_id()
            self.children.append(rel)      
    
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
    
    def __repr__(self):
        return f"Group(id={self.id}, name={self.name}, description={self.description})" 

class GroupRelation(Model):
    """
    Represents a relationship between two groups.

    A group relation defines a parent-child connection between two groups and
    indicates whether the child has admin privileges within the parent.

    Parameters
    ----------
    parent : Group
        The parent group in the relationship.
    child : Group
        The child group in the relationship.
    is_admin : bool, optional
        Whether the relationship includes admin privileges (default is False).

    Raises
    ------
    ValueError
        If either `parent` or `child` is not provided.
    """

    def __init__(
        self,
        parent: Optional[Group] = None,
        child: Optional[Group] = None,
        is_admin: bool = False,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.parent: Group = parent
        self.child: Group = child
        self.is_admin: bool = is_admin

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "GroupRelation":
        return cls(
            id=data.get("id"),
            parent=Group.from_dict(data["parent"]) if data["parent"] else None,
            child=Group.from_dict(data["child"]) if data["child"] else None,
            is_admin=data.get("isAdmin", False)
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

    def __repr__(self):
        return f"GroupRelation(parent={self.parent.__repr__()}, child={self.child.__repr__()}, is_admin={self.is_admin})"   
