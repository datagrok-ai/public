import pytest
import pandas as pd
from datagrok_api.models import DataConnection, DataSourceType


@pytest.fixture(scope='module')
def file_conn(grok):
    """Find a file-type connection on the server, skip if none exists."""
    conns = grok.connections.list()
    for c in conns:
        if DataSourceType.is_file_data_source(c.data_source):
            return c
    pytest.skip('No file-type connection available on the server')


def test_upload_and_download_csv(grok, file_conn, tmp_path):
    csv_file = tmp_path / 'test.csv'
    csv_file.write_text('col1,col2\n1,2\n3,4\n')

    grok.files.upload(file_conn, 'pytest_test.csv', str(csv_file))
    df = grok.files.download(file_conn, 'pytest_test.csv')
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ['col1', 'col2']
    assert len(df) == 2


def test_upload_and_download_binary(grok, file_conn, tmp_path):
    bin_file = tmp_path / 'test.bin'
    bin_file.write_bytes(b'\x00\x01\x02\x03')

    grok.files.upload(file_conn, 'pytest_test.bin', str(bin_file))
    content = grok.files.download(file_conn, 'pytest_test.bin')
    assert isinstance(content, bytes)
    assert content == b'\x00\x01\x02\x03'
