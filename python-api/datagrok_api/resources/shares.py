from typing import List, Union
from datagrok_api.http_client import HttpClient
from datagrok_api.models.group import Group
from datagrok_api.models.model import Model
from datagrok_api.models.share_response import ShareResponse


class SharesClient:
    """
    Client for managing shares of entities in the Datagrok platform.

    This client provides functionality to share entities (such as queries, data connections, scripts, etc.)
    with one or more groups or users, specifying access permissions.

    Parameters
    ----------
    client : HttpClient
        An instance of the HTTP client used to communicate with the Datagrok API.

    Methods
    -------
    share(id, groups, access='View') -> ShareResponse
        Shares the specified entity with given groups or users with the specified access level.

    Examples
    -------
    Prepare an entity to share

    >>> conn = DataConnection(
    ...     name="TestConn",
    ...     data_source=DatabaseDataSourceType.Postgres,
    ...     description="Test conn",
    ...     server="localhost",
    ...     port=5432,
    ...     db="datagrok",
    ...     credentials=Credentials(login="postgres", password="postgres")
    ... )
    >>> conn = grok.connections.save(conn, save_credentials=True)
    >>> query = grok.functions.create_query(connection=conn, query="select * from events", name="Events")

    Share the query with the DEVELOPERS group with View access
    >>> grok.shares.share(query, groups=[Group.DEVELOPERS()])
    """
    def __init__(self, client: HttpClient):
        self.client = client

    def share(self,
        id: Union[str, Model],
        groups: Union[str, List[Union[str, Group]]],
        access: str = "View") -> ShareResponse:
        '''
        Shares an entity.

        Parameters
        ----------
        id : str or Model
            Identifier of an entity: ID or Grok name, or a Model instance.
        groups : str or list of (str or Group)
            Comma-separated list of group/user names, or a list of names/Group instances.
        access : str
            Either 'View' or 'Edit'
        '''

        if isinstance(id, Model):
            id = id.id
        if isinstance(groups, list):
            group_ids = []
            for g in groups:
                if isinstance(g, Group):
                    group_ids.append(g.name)
                else:
                    group_ids.append(str(g))
            groups = ",".join(group_ids)

        id = id.replace(':', '.')
        endpoint = f"/public/v1/entities/{id}/shares"
        params = {
            'groups': groups,
            'access': access
        }
        response = self.client.post(endpoint, params=params)
        shareResponse = ShareResponse(response.json())
        shareResponse.raise_for_failure()
        return shareResponse
