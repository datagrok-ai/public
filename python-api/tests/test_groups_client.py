import uuid
from datagrok_api.models import Group


def _unique_name():
    return f'pytest_{uuid.uuid4().hex[:8]}'


def test_current_group(grok):
    group = grok.groups.current()
    assert isinstance(group, Group)
    assert group.id is not None


def test_list_groups(grok):
    groups = grok.groups.list()
    assert len(groups) > 0
    assert all(isinstance(g, Group) for g in groups)


def test_find_groups(grok):
    results = grok.groups.find('All users')
    assert len(results) >= 1
    assert any(g.name == 'All users' for g in results)


def test_save_and_delete_group(grok):
    name = _unique_name()
    group = Group(name=name, description='integration test')
    saved = grok.groups.save(group)
    try:
        assert saved.id is not None
        assert saved.name == name
    finally:
        grok.groups.delete(saved)


def test_get_group(grok):
    name = _unique_name()
    group = Group(name=name, description='integration test')
    saved = grok.groups.save(group)
    try:
        fetched = grok.groups.get(saved.id)
        assert fetched.id == saved.id
        assert fetched.name == name
    finally:
        grok.groups.delete(saved)


def test_add_member(grok):
    my_group = grok.groups.current()
    parent = Group(name=_unique_name(), description='parent group')
    parent = grok.groups.save(parent)
    try:
        result = grok.groups.add_member(parent, my_group)
        assert result.id == parent.id
    finally:
        grok.groups.delete(parent)


def test_get_members(grok):
    all_users = Group.ALL_USERS()
    members = grok.groups.get_members(all_users)
    assert isinstance(members, list)
    assert all(isinstance(g, Group) for g in members)
