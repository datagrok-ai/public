import json
import pytest
import requests_mock
import pandas as pd
from io import StringIO

from datagrok_api.http_client import HttpClient
from datagrok_api.resources import TablesClient

@pytest.fixture
def mock_client():
    with requests_mock.Mocker() as m:
        yield HttpClient("http://localhost:8888", "mock-token"), m

def test_download_table(mock_client):
    client, m = mock_client
    csv_content = "col1,col2\nval1,val2\nval3,val4"
    table_name = "Admin:.MyTestTable"
    m.get(f"http://localhost:8888/public/v1/tables/{table_name.replace(':', '.')}",
          text=csv_content,
          headers={'Content-Type': 'text/plain'})

    tables_client = TablesClient(client)
    df = tables_client.download(table_name)

    expected_df = pd.read_csv(StringIO(csv_content))
    pd.testing.assert_frame_equal(df, expected_df)

def test_upload_table(mock_client):
    client, m = mock_client
    table_name = "namespace.project.table"
    input_df = pd.DataFrame({"col1": ["val1", "val3"], "col2": ["val2", "val4"]})

    # Mock the POST request, returning some JSON response
    response_json = {"id": "table123", "name": table_name}
    m.post(f"http://localhost:8888/public/v1/tables/{table_name.replace(':', '.')}",
           json=response_json,
           headers={'Content-Type': 'application/json'})

    tables_client = TablesClient(client)
    result = tables_client.upload(table_name, input_df)

    # Verify the response JSON is as expected
    assert result == response_json

    # Verify the posted data is CSV matching input_df
    last_request = m.last_request
    posted_csv = last_request.body.decode('utf-8')
    expected_csv = input_df.to_csv(index=False)
    assert posted_csv == expected_csv

def test_create_dashboard_without_layout(mock_client):
    client, m = mock_client
    dashboard_name = "namespace.dashboard"
    table_ids = "table1,table2"

    response_json = {"dashboard_id": "dash123", "name": dashboard_name}
    endpoint = f"http://localhost:8888/public/v1/dashboards/{dashboard_name.replace(':', '.')}/{table_ids.replace(':', '.')}"
    m.post(endpoint, json=response_json)

    tables_client = TablesClient(client)
    result = tables_client.create_dashboard(dashboard_name, table_ids)

    assert result == response_json
    last_request = m.last_request

def test_create_dashboard_with_layout(tmp_path, mock_client):
    client, m = mock_client
    dashboard_name = "namespace.dashboard"
    table_ids = "table1,table2"
    layout = {"widgets": [{"type": "chart", "id": "widget1"}]}

    # Write layout JSON to temp file
    layout_file = tmp_path / "layout.json"
    layout_file.write_text(json.dumps(layout))

    response_json = {"dashboard_id": "dash123", "name": dashboard_name}
    endpoint = f"http://localhost:8888/public/v1/dashboards/{dashboard_name.replace(':', '.')}/{table_ids.replace(':', '.')}"
    m.post(endpoint, json=response_json)

    tables_client = TablesClient(client)
    result = tables_client.create_dashboard(dashboard_name, table_ids, layout_filename=str(layout_file))

    assert result == response_json
    last_request = m.last_request
    assert last_request.json() == layout
