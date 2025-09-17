import pytest
import requests_mock
from datagrok_api.http_client import HttpClient
from datagrok_api.models import DataConnection, DatabaseDataSourceType
from datagrok_api.models import ScriptLanguage
from datagrok_api.resources import SharesClient
from datagrok_api.models import Func
from datagrok_api.models import DataQuery
from datagrok_api.models import Script
from datagrok_api.models import Group
from datagrok_api.models import ShareResponse

@pytest.fixture
def mock_client():
    with requests_mock.Mocker() as m:
        yield HttpClient("http://localhost:8888", "mock-token"), m

def sample_func():
    return Func(id="func123", name="MyFunc", namespace="Namespace:Project:")

def sample_query():
    return DataQuery(id="query123", name="MyQuery", query="select 1", namespace="Namespace:Project:", connection=DataConnection(name="TestConn", data_source=DatabaseDataSourceType.Postgres))

def sample_script():
    return Script(id="script123", name="MyScript", script="#name:MyScript\nprint('Hello World')", namespace="Namespace:Project:", language=ScriptLanguage.Python)

def sample_group(name="group1"):
    return Group(id="group123", name=name)

def test_share_with_string_id_and_string_group(mock_client):
    client, m = mock_client
    shares_client = SharesClient(client)
    m.post("http://localhost:8888/public/v1/entities/entity123/shares", json={"success": True})

    resp = shares_client.share("entity123", "groupA", access="View")
    assert isinstance(resp, ShareResponse)

def test_share_with_model_id_and_string_group(mock_client):
    client, m = mock_client
    shares_client = SharesClient(client)
    func = sample_func()
    m.post("http://localhost:8888/public/v1/entities/func123/shares", json={"success": True})

    resp = shares_client.share(func, "groupA", access="Edit")
    assert isinstance(resp, ShareResponse)

def test_share_with_string_id_and_list_of_strings(mock_client):
    client, m = mock_client
    shares_client = SharesClient(client)
    m.post("http://localhost:8888/public/v1/entities/entity123/shares", json={"success": True})

    resp = shares_client.share("entity123", ["groupA", "groupB"], access="View")
    assert isinstance(resp, ShareResponse)

def test_share_with_string_id_and_list_of_groups(mock_client):
    client, m = mock_client
    shares_client = SharesClient(client)
    group1 = sample_group("groupA")
    group2 = sample_group("groupB")
    m.post("http://localhost:8888/public/v1/entities/entity123/shares", json={"success": True})

    resp = shares_client.share("entity123", [group1, group2], access="Edit")
    assert isinstance(resp, ShareResponse)

def test_share_with_query_and_mixed_groups(mock_client):
    client, m = mock_client
    shares_client = SharesClient(client)
    query = sample_query()
    group1 = sample_group("groupA")
    m.post("http://localhost:8888/public/v1/entities/query123/shares", json={"success": True})

    resp = shares_client.share(query, ["groupB", group1], access="View")
    assert isinstance(resp, ShareResponse)

def test_share_with_script_and_empty_groups(mock_client):
    client, m = mock_client
    shares_client = SharesClient(client)
    script = sample_script()
    m.post("http://localhost:8888/public/v1/entities/script123/shares", json={"success": True})

    resp = shares_client.share(script, "", access="View")
    assert isinstance(resp, ShareResponse)
