"""Client output schemas. The provider constrains the reply to what its grammar supports; the reply is then checked
here against the whole schema, in the session's own process, so a pattern that backtracks stalls only this session,
which the runtime's idle watchdog ends."""
import json
from typing import Any

from jsonschema import Draft7Validator
from jsonschema.exceptions import SchemaError
from jsonschema.protocols import Validator
from jsonschema.validators import validator_for

from .protocol import JsonObject


def _validator_class(schema: JsonObject) -> type[Validator]:
    # Draft-07 when the schema declares no `$schema`: the draft most client schemas are written in.
    validator: type[Validator] = validator_for(schema, default=Draft7Validator)
    return validator


def schema_error(schema: Any) -> str | None:
    if not isinstance(schema, dict):
        return 'it is not a JSON object'
    try:
        _validator_class(schema).check_schema(schema)
    except SchemaError as e:
        return e.message
    return None


def check(schema: JsonObject, text: str) -> tuple[Any, str | None]:
    """(value, None) for a valid reply, else (None, what is wrong)."""
    try:
        value = json.loads(text)
    except ValueError:
        return None, 'the reply is not a JSON value'
    errors = sorted(_validator_class(schema)(schema).iter_errors(value), key=lambda e: list(e.path))
    if errors:
        return None, '; '.join(f"{'/'.join(map(str, e.path)) or 'reply'}: {e.message}" for e in errors[:5])
    return value, None
