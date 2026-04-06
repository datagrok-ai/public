import uuid
from datagrok_api.models import DataConnection, DatabaseDataSourceType, Credentials, ShareResponse


def test_share_entity(grok):
    name = f'pytest_{uuid.uuid4().hex[:8]}'
    conn = DataConnection(
        name=name,
        data_source=DatabaseDataSourceType.Postgres,
        server='localhost',
        port=5432,
        db='test_db',
        credentials=Credentials(login='user', password='pass'),
    )
    conn = grok.connections.save(conn)
    try:
        resp = grok.shares.share(conn, 'All users', access='View')
        assert isinstance(resp, ShareResponse)
    finally:
        grok.connections.delete(conn)
