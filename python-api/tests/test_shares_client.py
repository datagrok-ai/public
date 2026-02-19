import uuid
from datagrok_api.models import Group, ShareResponse


def test_share_entity(grok):
    name = f'pytest_{uuid.uuid4().hex[:8]}'
    group = Group(name=name, description='share test')
    group = grok.groups.save(group)
    try:
        resp = grok.shares.share(group, 'All users', access='View')
        assert isinstance(resp, ShareResponse)
    finally:
        grok.groups.delete(group)
