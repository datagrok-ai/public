from datetime import datetime
from typing import Optional
from datagrok_api.models.data_connection import DataConnection
from datagrok_api.models.group import Group
from datagrok_api.models.model import NamedModel
from enum import Enum


class UserStatus(str, Enum):
    new = "new"
    active = "active"
    blocked = "blocked"
    guest = "guest"


class User(NamedModel):
    """
    Represents a user in Datagrok with personal and account-related details.

    Parameters
    ----------
    first_name : Optional[str]
        User's first name.
    last_name : Optional[str]
        User's last name.
    login : Optional[str]
        Username or login identifier.
    email : Optional[str]
        User's email address.
    status : Optional[UserStatus]
        Status of the user account (e.g., active, inactive). Defaults to UserStatus.new.
    storage : Optional[DataConnection]
        Data connection representing user's storage.
    group : Optional[Group]
        Group associated with the user.
    **kwargs
        Additional keyword arguments passed to the base NamedModel constructor.

    Examples
    --------
    Create a new user with login and email

    >>> user = User(login="jdoe", email="jdoe@example.com", first_name="John", last_name="Doe")
    """
    def __init__(self, 
                 first_name: Optional[str] = None, 
                 last_name: Optional[str] = None, 
                 login: Optional[str] = None, 
                 email: Optional[str] = None, 
                 status: Optional[UserStatus] = None, 
                 storage: Optional[DataConnection] = None, 
                 group: Optional[Group] =None,
                 **kwargs):
        super().__init__(**kwargs)
        self.first_name = first_name
        self.last_name = last_name
        self.login = login
        self.email = email
        self.status = status or UserStatus.new
        self.storage = storage
        self.group = group


    def to_dict(self):
        return {
            "id": self.id,
            "name": self.name,
            "friendlyName": self.friendly_name,
            "firstName": self.first_name,
            "lastName": self.last_name,
            "email": self.email,
            "login": self.login,
            "status": self.status.value,
            "storage": None if self.storage is None else self.storage.to_dict(),
            "group": None if self.group is None else self.group.to_dict(),
            "createdOn": self.created_on.isoformat() if self.created_on else None,
            "updatedOn": self.updated_on.isoformat() if self.updated_on else None,
            "namespace": self.namespace
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> "User":
        return cls(
            id=data.get("id"),
            name=data.get("name"),
            friendly_name=data.get("friendlyName"),
            first_name=data.get("firstName"),
            last_name=data.get("lastName"),
            login=data.get("login"),
            email=data.get("email"),
            status=UserStatus(data.get("status")),
            storage=DataConnection.from_dict(data["storage"]) if data.get("storage") else None,
            group=Group.from_dict(data["group"]) if data.get("group") else None,
            created_on=datetime.fromisoformat(data["createdOn"]) if data.get("createdOn") else None,
            updated_on=datetime.fromisoformat(data["updatedOn"]) if data.get("updatedOn") else None,
            namespace=data.get("namespace")
        )
