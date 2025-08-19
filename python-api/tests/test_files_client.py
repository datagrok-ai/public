import pytest
import pandas as pd
import requests_mock

from datagrok_api.http_client import HttpClient
from datagrok_api.models import DataConnection, FileDataSourceType
from datagrok_api.resources import FilesClient


@pytest.fixture
def mock_client():
    with requests_mock.Mocker() as m:
        yield HttpClient("http://localhost:8888", "mock-token"), m


def sample_connection():
    return DataConnection(
        id="test-id",
        name="TestConn",
        namespace="Admin:",
        data_source=FileDataSourceType.S3,
        bucket='datagrok'
    )


def test_download_csv(mock_client):
    client, m = mock_client
    files_client = FilesClient(client)

    csv_content = "col1,col2\n1,2\n3,4\n"
    endpoint = "http://localhost:8888/public/v1/files/Admin.TestConn/test.csv"
    m.get(endpoint, content=csv_content.encode('utf-8'))

    df = files_client.download("Admin:TestConn", "test.csv")
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["col1", "col2"]
    assert df.iloc[0, 0] == 1
    assert df.iloc[1, 1] == 4

    df = files_client.download(sample_connection(), "test.csv")
    assert isinstance(df, pd.DataFrame)


def test_download_binary(mock_client):
    client, m = mock_client
    files_client = FilesClient(client)

    binary_content = b'\x00\x01\x02\x03'
    endpoint = "http://localhost:8888/public/v1/files/Admin.TestConn/test.bin"
    m.get(endpoint, content=binary_content)

    content = files_client.download(sample_connection(), "test.bin")
    assert isinstance(content, bytes)
    assert content == binary_content


def test_upload_file(tmp_path, mock_client):
    client, m = mock_client
    files_client = FilesClient(client)

    # Prepare a small test file
    test_file = tmp_path / "upload.txt"
    test_file.write_text("hello world")

    endpoint = "http://localhost:8888/public/v1/files/Admin.TestConn/upload.txt"

    # Mock POST endpoint - here we don't care about the response content, just that it was called
    m.post(endpoint, status_code=200)

    files_client.upload("Admin:TestConn", "upload.txt", str(test_file))

    # Assert the POST was called exactly once with correct headers and data
    assert m.called
    req = m.request_history[0]
    assert req.method == "POST"
    assert req.url == endpoint
    assert req.headers.get('Content-Type') == 'application/octet-stream'
