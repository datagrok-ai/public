import pandas as pd


def test_upload_and_download(grok):
    df = pd.DataFrame({'name': ['Alice', 'Bob'], 'age': [30, 25]})
    result = grok.tables.upload('pytest_test_table', df)
    assert isinstance(result, dict)
    table_id = result.get('ID') or result.get('id')
    assert table_id is not None

    downloaded = grok.tables.download(table_id)
    assert isinstance(downloaded, pd.DataFrame)
    assert list(downloaded.columns) == ['name', 'age']
    assert len(downloaded) == 2
