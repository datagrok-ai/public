import uuid
from datagrok_api.models import Func, Script


def _unique_name():
    return f'pytest_{uuid.uuid4().hex[:8]}'


def test_list_functions(grok):
    funcs = grok.functions.list()
    assert isinstance(funcs, list)
    assert all(isinstance(f, Func) for f in funcs)


def test_list_scripts(grok):
    scripts = grok.functions.list_scripts()
    assert isinstance(scripts, list)
    assert all(isinstance(s, Script) for s in scripts)


def test_create_and_get_script(grok):
    name = _unique_name()
    script = grok.functions.create_script(script="print('hello')", name=name)
    assert isinstance(script, Script)
    assert script.id is not None

    fetched = grok.functions.get(script.id)
    assert fetched.id == script.id
    assert fetched.name == name
