from typing import Union
from urllib.parse import quote
from datagrok_api.http_client import HttpClient
from datagrok_api.models.data_connection import Credentials


class CredentialsClient:
    """
    Client for managing credentials for any Datagrok entity.

    Credentials can be associated with data connections, packages,
    or other entity types.

    Parameters
    ----------
    client : HttpClient
        An initialized HTTP client used for making API requests.

    Examples
    --------
    Get credentials for a connection

    >>> creds = grok.credentials.for_entity(connection.id)
    >>> print(creds.parameters)

    Save credentials for a package

    >>> creds = grok.credentials.for_entity(package.id)
    >>> creds.parameters['apiKey'] = 'new-key'
    >>> grok.credentials.save(creds)
    """
    def __init__(self, client: HttpClient):
        self.client = client

    def for_entity(self, entity_id: str) -> Credentials:
        """Get credentials for an entity by its ID.

        Parameters
        ----------
        entity_id : str
            The entity ID (UUID) to retrieve credentials for.

        Returns
        -------
        Credentials
            The credentials associated with the entity.
        """
        endpoint = f"/public/v1/credentials/for/{quote(entity_id)}"
        response = self.client.get(endpoint)
        return Credentials.from_dict(response.json())

    def save(self, credentials: Credentials) -> None:
        """Save or update credentials.

        The credentials must have ``entity_bind_id`` set to associate
        them with an entity. Typically, first call ``for_entity()`` to
        load existing credentials, modify parameters, then save.

        Parameters
        ----------
        credentials : Credentials
            The credentials object to save.
        """
        endpoint = "/public/v1/credentials"
        self.client.post(endpoint, json=credentials.to_dict())

    def delete(self, credentials: Union[Credentials, str]) -> None:
        """Delete credentials by object or ID.

        Parameters
        ----------
        credentials : Credentials or str
            The credentials object or credentials ID to delete.
        """
        cred_id = credentials.id if isinstance(credentials, Credentials) else credentials
        endpoint = f"/public/v1/credentials/{quote(cred_id)}"
        self.client.delete(endpoint)
