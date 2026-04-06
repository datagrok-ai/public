import uuid
from datagrok_api.models import DataConnection, DatabaseDataSourceType, Credentials


def _unique_name():
    return f'pytest_{uuid.uuid4().hex[:8]}'


def test_credentials_for_connection(grok):
    name = _unique_name()
    conn = DataConnection(
        name=name,
        data_source=DatabaseDataSourceType.Postgres,
        server='localhost',
        port=5432,
        db='test_db',
        credentials=Credentials(login='test_user', password='test_pass'),
    )
    saved = grok.connections.save(conn, save_credentials=True)
    try:
        creds = grok.credentials.for_entity(saved.id)
        assert creds is not None
        assert creds.entity_bind_id is not None
        assert creds.group is not None
    finally:
        grok.connections.delete(saved)


def test_credentials_save(grok):
    name = _unique_name()
    conn = DataConnection(
        name=name,
        data_source=DatabaseDataSourceType.Postgres,
        server='localhost',
        port=5432,
        db='test_db',
        credentials=Credentials(login='test_user', password='test_pass'),
    )
    saved = grok.connections.save(conn, save_credentials=True)
    try:
        creds = grok.credentials.for_entity(saved.id)
        creds.parameters['login'] = 'updated_user'
        grok.credentials.save(creds)

        updated = grok.credentials.for_entity(saved.id)
        assert updated is not None
    finally:
        grok.connections.delete(saved)


def test_package_credentials_roundtrip(grok):
    """Get ApiTests package, set credentials, verify, then delete."""
    pkg = grok.packages.get('ApiTests')
    assert pkg is not None

    # Clean up any pre-existing credentials
    try:
        old_creds = grok.credentials.for_entity(pkg.id)
        grok.credentials.delete(old_creds)
    except Exception:
        pass

    token_value = f'test_token_{uuid.uuid4().hex[:8]}'

    creds = Credentials(entity_bind_id=pkg.id, external_token=token_value)
    grok.credentials.save(creds)

    # Retrieve and verify the value was stored
    saved_creds = grok.credentials.for_entity(pkg.id)
    assert saved_creds.parameters.get('external_token') == token_value
    assert saved_creds.group is not None

    # Clean up
    grok.credentials.delete(saved_creds)
