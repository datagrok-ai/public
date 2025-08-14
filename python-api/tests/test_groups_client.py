import pytest
import requests_mock
from datagrok_api.http_client import HttpClient
from datagrok_api.models import Group
from datagrok_api.models import User
from datagrok_api.resources import GroupsClient

@pytest.fixture
def mock_client():
    with requests_mock.Mocker() as m:
        yield HttpClient("http://localhost:8888", "mock-token"), m

def sample_group():
    return Group(id="group1", name="TestGroup", description="Testing", friendly_name="Test Group")

def sample_user():
    return User(id="user1", first_name="John", last_name="Doe", login="johnydoe", email="johnydoe@example.com")

def test_find_groups(mock_client):
    client, m = mock_client
    groups = [sample_group().to_dict()]
    m.get("http://localhost:8888/public/v1/groups/lookup", json=groups)
    result = GroupsClient(client).find("TestGroup")
    assert len(result) == 1
    assert result[0].id == "group1"

def test_get_group(mock_client):
    client, m = mock_client
    group = sample_group().to_dict()
    m.get("http://localhost:8888/public/v1/groups/group1", json=group)
    result = GroupsClient(client).get(group["id"])
    assert result.id == "group1"
    assert result.friendly_name == "Test Group"

def test_save_group(mock_client):
    client, m = mock_client
    group = sample_group()
    saved_group = group.to_dict()
    m.post("http://localhost:8888/public/v1/groups", json=saved_group)
    result = GroupsClient(client).save(group)
    assert result.id == group.id
    assert result.name == group.name

def test_delete_group_by_id(mock_client):
    client, m = mock_client
    group_id = "group1"
    m.delete(f"http://localhost:8888/public/v1/groups/{group_id}")
    GroupsClient(client).delete(group_id)
    # no exception means success

def test_delete_group_by_object(mock_client):
    client, m = mock_client
    group = sample_group()
    m.delete(f"http://localhost:8888/public/v1/groups/{group.id}")
    GroupsClient(client).delete(group)
    # no exception means success

def test_list_groups_with_options(mock_client):
    client, m = mock_client
    groups = [sample_group().to_dict()]
    url = "http://localhost:8888/public/v1/groups"
    m.get(url, json=groups)
    result = GroupsClient(client).list(smart_filter="name='Test'", include_personal=True, include_members=True, include_memberships=True)
    assert len(result) == 1
    assert result[0].id == "group1"

def test_add_member_with_group_and_user(mock_client):
    client, m = mock_client
    groups_client = GroupsClient(client)

    parent_group = sample_group()
    child_group = Group(id="group2", name="Child Group")
    user = sample_user()

    # Mock find to return exactly one group for parent and child
    m.get("http://localhost:8888/public/v1/groups/lookup", json=[parent_group.to_dict()])
    m.get("http://localhost:8888/public/v1/groups/lookup", json=[child_group.to_dict()])
    # Mock find for user name to find personal group
    personal_group = Group(id="personal1", name=user.name, personal=True)
    m.get(f"http://localhost:8888/public/v1/groups/lookup", json=[personal_group.to_dict()])

    # Mock save call
    m.post("http://localhost:8888/public/v1/groups", json=parent_group.to_dict())

    # Patch add_member method to avoid implementation detail and just add child
    def dummy_add_member(self, child, is_admin=False):
        assert child.id == personal_group.id or child.id == child_group.id

    parent_group.add_member = dummy_add_member.__get__(parent_group, Group)

    # Test with child as User (should resolve to personal group)
    result = groups_client.add_member(parent_group, user)
    assert result.id == parent_group.id

def test_get_members_and_memberships(mock_client):
    client, m = mock_client
    groups_client = GroupsClient(client)
    group = sample_group()

    members = [sample_group().to_dict()]
    memberships = [sample_group().to_dict()]

    m.get(f"http://localhost:8888/public/v1/groups/{group.id}/members", json=members)
    m.get(f"http://localhost:8888/public/v1/groups/{group.id}/memberships", json=memberships)

    result_members = groups_client.get_members(group)
    result_members_admin = groups_client.get_members(group, admin=True)
    result_memberships = groups_client.get_memberships(group)
    result_memberships_admin = groups_client.get_memberships(group, admin=True)

    assert all(isinstance(g, Group) for g in result_members)
    assert all(isinstance(g, Group) for g in result_members_admin)
    assert all(isinstance(g, Group) for g in result_memberships)
    assert all(isinstance(g, Group) for g in result_memberships_admin)

def test_current_group(mock_client):
    client, m = mock_client
    current_group = sample_group().to_dict()
    m.get("http://localhost:8888/public/v1/groups/current", json=current_group)
    result = GroupsClient(client).current()
    assert result.id == current_group["id"]
