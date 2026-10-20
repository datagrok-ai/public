"""Hermes takes hooks only from plugins: this one registers the hooks of the session running in this process."""
from grokky_agent import session
from grokky_agent.hermes import PluginContext


def register(ctx: PluginContext) -> None:
    assert session.current is not None, 'the session is created before Hermes loads its plugins'
    session.current.register_hooks(ctx)
