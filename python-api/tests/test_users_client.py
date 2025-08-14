import pytest
import requests_mock
from datagrok_api.http_client import HttpClient
from datagrok_api.resources import UsersClient  # adjust import to your project structure
from datagrok_api.models import User

@pytest.fixture
def mock_client():
    with requests_mock.Mocker() as m:
        yield HttpClient("http://localhost:8888", "mock-token"), m

def sample_user():
    return User(
        id="user123",
        name="John Doe",
        first_name="John",
        last_name="Doe",
        login="johndoe",
        email="johndoe@example.com"
    )

def test_get_user(mock_client):
    client, m = mock_client
    user = sample_user()
    user_id = user.id
    m.get(f"http://localhost:8888/public/v1/users/{user_id}", json=user.to_dict())

    users_client = UsersClient(client)
    result = users_client.get(user_id)
    assert result.id == user.id
    assert result.login == user.login
    assert result.email == user.email

def test_save_user(mock_client):
    client, m = mock_client
    user = sample_user()
    m.post("http://localhost:8888/public/v1/users", json=user.to_dict())

    users_client = UsersClient(client)
    result = users_client.save(user)
    assert result.id == user.id
    assert result.login == user.login

def test_current_user(mock_client):
    client, m = mock_client
    user = sample_user()
    m.get("http://localhost:8888/public/v1/users/current", json=user.to_dict())

    users_client = UsersClient(client)
    result = users_client.current()
    assert result.id == user.id

def test_list_users(mock_client):
    client, m = mock_client
    user1 = sample_user()
    user2 = User(id="user456", name="Jane Smith", login="janesmith", email="jane@example.com")
    m.get("http://localhost:8888/public/v1/users", json=[user1.to_dict(), user2.to_dict()])

    users_client = UsersClient(client)
    results = users_client.list()
    assert len(results) == 2
    assert any(u.id == user1.id for u in results)
    assert any(u.id == user2.id for u in results)

def test_block_user(mock_client):
    client, m = mock_client
    user = sample_user()
    m.post("http://localhost:8888/public/v1/users/block", json={})

    users_client = UsersClient(client)
    users_client.block(user)
    # Optionally check last request body if you want
    last_request = m.last_request
    assert last_request.json() == user.to_dict()

def test_unblock_user(mock_client):
    client, m = mock_client
    user = sample_user()
    m.post("http://localhost:8888/public/v1/users/unblock", json={})

    users_client = UsersClient(client)
    users_client.unblock(user)
    last_request = m.last_request
    assert last_request.json() == user.to_dict()
