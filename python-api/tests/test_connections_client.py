import pytest
import requests_mock

from datagrok_api.http_client import HttpClient
from datagrok_api.models import DataConnection, DatabaseDataSourceType, Credentials
from datagrok_api.resources import ConnectionsClient

@pytest.fixture
def mock_client():
    with requests_mock.Mocker() as m:
        yield HttpClient("http://localhost:8888", "mock-token"), m

def sample_connection():
    creds = Credentials(login="datagrok", password="dg1234")
    return DataConnection(
        id="test-id",
        name="Test Connection",
        data_source=DatabaseDataSourceType.Postgres,
        server="db.datagrok.ai",
        port=4332,
        db="test",
        credentials=creds
    )

def test_save_connection(mock_client):
    client, m = mock_client
    conn = sample_connection()
    m.post("http://localhost:8888/public/v1/connections", json=conn.to_dict())
    result = ConnectionsClient(client).save(conn)
    assert result.name == conn.name
    assert result.data_source == conn.data_source

def test_delete_connection(mock_client):
    client, m = mock_client
    conn = sample_connection()
    m.delete(f"http://localhost:8888/public/v1/connections/{conn.id.replace(':', '.')}", status_code=204)
    ConnectionsClient(client).delete(conn) 

def test_list_connections(mock_client):
    client, m = mock_client
    conn = sample_connection()
    m.get("http://localhost:8888/public/v1/connections", json=[conn.to_dict()])
    conns = ConnectionsClient(client).list()
    assert len(conns) == 1
    assert conns[0].name == conn.name

def test_test_connection_success(mock_client):
    client, m = mock_client
    conn = sample_connection()
    m.post("http://localhost:8888/public/v1/connections/test", json="ok")
    ConnectionsClient(client).test(conn)

