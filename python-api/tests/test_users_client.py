from datagrok_api.models import User


def test_current_user(grok):
    user = grok.users.current()
    assert isinstance(user, User)
    assert user.id is not None
    assert user.login is not None


def test_get_user(grok):
    me = grok.users.current()
    user = grok.users.get(me.id)
    assert user.id == me.id
    assert user.login == me.login


def test_list_users(grok):
    users = grok.users.list()
    assert len(users) > 0
    assert all(isinstance(u, User) for u in users)
    me = grok.users.current()
    assert any(u.id == me.id for u in users)
