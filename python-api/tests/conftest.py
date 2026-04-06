import os
import pytest
import requests
from datagrok_api import DatagrokClient


def pytest_addoption(parser):
    parser.addoption("--url", default=os.environ.get("DATAGROK_API_URL", "http://localhost:8082"))
    parser.addoption("--key", default=os.environ.get("DATAGROK_API_KEY", "unit_test_token"))


@pytest.fixture(scope="session")
def grok(request):
    url = request.config.getoption("--url")
    key = request.config.getoption("--key")
    try:
        requests.get(f"{url}/info/server", timeout=3)
    except requests.ConnectionError:
        pytest.skip(f"Datagrok server not reachable at {url}")
    client = DatagrokClient(base_url=url, api_key=key)
    yield client
    client.close()
