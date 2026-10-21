from abc import ABC, abstractmethod

from ..hermes import Connection
from ..protocol import JsonObject


class ProviderError(ValueError):
    """The container credentials cannot be used with this provider."""


class Provider(ABC):
    """A model API the agent runs on: how Hermes reaches it, and how a turn's system prompt and output format go on
    a request."""

    @abstractmethod
    def connection(self) -> Connection: ...

    @abstractmethod
    def patch_request(self, request: JsonObject, system_prompt: str, output_schema: JsonObject | None) -> JsonObject:
        """The outgoing request with the turn's system prompt (none when empty) and output format."""

    def refused_request(self) -> JsonObject:
        """Sent instead of a request the session cannot vouch for: the API rejects it, and nothing of the
        conversation leaves the process."""
        return {'model': 'grokky-refused', 'max_tokens': 1, 'messages': [{'role': 'user', 'content': 'refused'}]}

    def error_hint(self, status: int | None) -> str:
        """Added to a failed call's message, e.g. what to check for a rejected key."""
        return ''
