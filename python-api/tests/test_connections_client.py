import uuid
from datagrok_api.models import DataConnection, DatabaseDataSourceType, Credentials


def _unique_name():
    return f'pytest_{uuid.uuid4().hex[:8]}'


def test_list_connections(grok):
    conns = grok.connections.list()
    assert isinstance(conns, list)
    assert all(isinstance(c, DataConnection) for c in conns)


def test_save_get_delete_connection(grok):
    name = _unique_name()
    conn = DataConnection(
        name=name,
        data_source=DatabaseDataSourceType.Postgres,
        server='localhost',
        port=5432,
        db='test_db',
        credentials=Credentials(login='user', password='pass'),
    )
    saved = grok.connections.save(conn)
    try:
        assert saved.id is not None
        assert saved.name == name

        fetched = grok.connections.get(saved.id)
        assert fetched.id == saved.id
        assert fetched.name == name
    finally:
        grok.connections.delete(saved)
