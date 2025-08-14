import pytest
import requests_mock

from datagrok_api.http_client import HttpClient
from datagrok_api.resources import FunctionsClient
from datagrok_api.models import Func, DataQuery, Script, ScriptLanguage, Credentials
from datagrok_api.models import DataConnection, DatabaseDataSourceType


@pytest.fixture
def mock_client():
    with requests_mock.Mocker() as m:
        yield HttpClient("http://localhost:8888", "mock-token"), m

def sample_func():
    return Func(id="123", namespace="Admin:", name="TestFunc")

def sample_query():
    return DataQuery(id="q1", name="Query1", query="SELECT 1", connection=sample_connection(), namespace="Admin:")

def sample_script():
    return Script(id="s1", name="Script1", script="print('hi')", language=ScriptLanguage.Python, namespace="Admin:")

def sample_connection():
    creds = Credentials(login="datagrok", password="dg1234")
    return DataConnection(
        id="test-id",
        name="TestConnection",
        friendly_name="Test Connection",
        data_source=DatabaseDataSourceType.Postgres,
        server="db.datagrok.ai",
        port=4332,
        db="test",
        credentials=creds,
        namespace="Admin:"
    )

def test_call_function(mock_client):
    client, m = mock_client
    func_name = "package:func"
    parameters = {"x": 1}
    response_data = {"result": 2}
    endpoint = "http://localhost:8888/public/v1/functions/package.func/call"
    
    m.post(endpoint, json=response_data)

    fc = FunctionsClient(client)
    result = fc.call(func_name, parameters)
    assert result == response_data

def test_call_query(mock_client):
    client, m = mock_client
    query = sample_query()
    parameters = {"x": 1}
    response_data = {"result": 2}
    endpoint = "http://localhost:8888/public/v1/functions/Admin.Query1/call"
    
    m.post(endpoint, json=response_data)

    fc = FunctionsClient(client)
    result = fc.call(query, parameters)
    assert result == response_data

def test_call_script(mock_client):
    client, m = mock_client
    script = sample_script()
    parameters = {"x": 1}
    response_data = {"result": 2}
    endpoint = "http://localhost:8888/public/v1/functions/Admin.Script1/call"
    
    m.post(endpoint, json=response_data)

    fc = FunctionsClient(client)
    result = fc.call(script, parameters)
    assert result == response_data

def test_get_function(mock_client):
    client, m = mock_client
    func = sample_func()
    endpoint = f"http://localhost:8888/public/v1/functions/Admin.TestFunc"
    
    m.get(endpoint, json=func.to_dict())

    fc = FunctionsClient(client)
    fetched = fc.get("Admin:TestFunc")
    assert fetched.name == func.name
    assert fetched.id == func.id

def test_list_functions(mock_client):
    client, m = mock_client
    func = sample_func()
    endpoint = "http://localhost:8888/public/v1/functions"
    
    m.get(endpoint, json=[func.to_dict()])

    fc = FunctionsClient(client)
    funcs = fc.list()
    assert len(funcs) == 1
    assert funcs[0].id == func.id   


def test_list_queries(mock_client):
    client, m = mock_client
    query = sample_query()
    endpoint = "http://localhost:8888/public/v1/functions"
    
    m.get(endpoint, json=[query.to_dict()])

    fc = FunctionsClient(client)
    queries = fc.list_queries()
    assert len(queries) == 1
    assert queries[0].name == query.name

def test_list_scripts(mock_client):
    client, m = mock_client
    script = sample_script()
    endpoint = "http://localhost:8888/public/v1/functions"
    
    m.get(endpoint, json=[script.to_dict()])

    fc = FunctionsClient(client)
    scripts = fc.list_scripts()
    assert len(scripts) == 1
    assert scripts[0].name == script.name

def test_create_query(mock_client):
    client, m = mock_client
    query = sample_query()
    endpoint = "http://localhost:8888/public/v1/functions"

    m.post(endpoint, json=query.to_dict())

    fc = FunctionsClient(client)
    result = fc.create_query(query.connection, query.query, query.name)
    assert isinstance(result, DataQuery)
    assert result.name == query.name

def test_create_script(mock_client):
    client, m = mock_client
    script = sample_script()
    endpoint = "http://localhost:8888/public/v1/functions"

    m.post(endpoint, json=script.to_dict())

    fc = FunctionsClient(client)
    result = fc.create_script(script.script, script.name)
    assert isinstance(result, Script)
    assert result.script == script.script     
    