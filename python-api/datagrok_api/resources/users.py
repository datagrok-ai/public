from typing import List, Optional, Union
from datagrok_api.http_client import HttpClient
from datagrok_api.models.user import User
from urllib.parse import quote

class UsersClient:
    """
    Client for managing users in the Datagrok platform.

    Provides methods to retrieve, update, list, block, and unblock users,
    as well as to get information about the current authenticated user.

    Parameters
    ----------
    client : HttpClient
        An instance of the HTTP client used to communicate with the Datagrok API.

     Examples
    --------
    Get user by login

    >>> from datagrok_api import DatagrokClient
    >>> from datagrok_api.models import Func, Script, ScriptLanguage, DataQuery
    >>> grok = DatagrokClient(base_url="https://public.datagrok.ai/api", api_key="Bearer <your-token>")
    >>> john_doe = grok.users.get("johndoe")
    >>> print(john_doe.first_name)
    John

    Get current authenticated user
    >>> me = grok.users.current()
    """
    def __init__(self, client: HttpClient):
        self.client = client

    def get(self, id: str) -> User:
        """Get detailed information about a specific user.
        
        Parameters
        ----------
        id : str
            ID or login of the user to retrieve
            
        Returns
        -------
        User
            User instance with detailed information
        """
        endpoint = f"/public/v1/users/{quote(id)}"
        response = self.client.get(endpoint)
        return User.from_dict(response.json())

    def save(self, user: User) -> User:
        """Save or update User.
        
        Parameters
        ----------
        user : User
            The user to save or update
            
        Returns
        -------
        User
            The saved user instance
        """
        user.ensure_id() 
        endpoint = "/public/v1/users"
        response = self.client.post(endpoint, json=user.to_dict())
        return User.from_dict(response.json())

    def current(self) -> User:
        """Returns the current authenticated user.
        
        Returns
        -------
        User
            User instance representing the current user
        """
        endpoint = f"/public/v1/users/current"
        response = self.client.get(endpoint)
        return User.from_dict(response.json())
    
    def list(self, smart_filter: Optional[str] = None, include_group: bool = False) -> List[User]:
        """List users from Datagrok with optional filtering and inclusion of related data.
        
        Parameters
        ----------
        smart_filter : Optional[str]
            Optional smart search filter to apply
        include_group : bool, default=False
            Whether to include personal group in the results 

        Returns
        -------
        List[User]
            List of User instances matching the criteria
        """
        
        params = {}
        if smart_filter:
            params["text"] = smart_filter

        include_parts = []
        if include_group:
            include_parts.append("group")

        if include_parts:
            params["include"] = ",".join(include_parts)

        endpoint = "/public/v1/users"
        response = self.client.get(endpoint, params=params)
        
        return [User.from_dict(data) for data in response.json()]
    
    def block(self, user: User):
        endpoint = "/public/v1/users/block"
        self.client.post(endpoint, json=user.to_dict())

    def unblock(self, user: User):
        endpoint = "/public/v1/users/unblock"
        self.client.post(endpoint, json=user.to_dict())
            