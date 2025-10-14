from datetime import datetime
from typing import Optional
from datagrok_api.models.data_connection import DataConnection
from datagrok_api.models.func import Func, FuncParam

from enum import Enum


class ScriptLanguage(str, Enum):
    """
    Enumeration of supported scripting languages for Datagrok `Script` functions.
    Currently only server-side functions can be called from Python API.

    A `Script` function can be written in various programming languages
    supported by the platform. This enum defines the set of valid language
    identifiers used for both serialization and execution.

    Attributes
    -------
    Grok : str
        Grok's native scripting language.
    Julia : str
        The Julia scientific computing language.
    Python : str
        Python programming language.
    R : str
        The R statistical computing language.
    NodeJs : str
        JavaScript runtime environment (Node.js).
    Octave : str
        GNU Octave, primarily for numerical computations.

    Examples
    --------
    >>> ScriptLanguage.Python.value
    ... 'python'
    >>> ScriptLanguage['Julia']
    ... <ScriptLanguage.Julia: 'julia'>
    """
    Grok = 'grok'
    Julia = 'julia'
    Python = 'python'
    R = 'r'
    NodeJs = 'nodejs'
    Octave = 'octave'

@Func.register_subclass("script")
class Script(Func):
    """
    Represents a script-based function in the Datagrok platform.

    `Script` is a specialized subclass of `Func` that encapsulates 
    executable code written in a supported scripting language.  
    Scripts can be run in various environments (server-side or client-side) 
    depending on the language and platform configuration.

    Attributes
    ----------
    script : str
        The source code of the script. Can contain special annotations that will 
        be parsed by platform into script paramerters.
    language : ScriptLanguage, optional
        The language in which the script is written. Defaults to `ScriptLanguage.Python`.
    **kwargs
        Additional arguments passed to the base `Func` class, such as:
            - id : Optional[str]
            - name : Optional[str]
            - friendly_name : Optional[str]
            - description : Optional[str]
            - created_on : Optional[datetime]
            - updated_on : Optional[datetime]
            - params : Optional[List[FuncParam]]
            - namespace : Optional[str]
            - tags : Optional[List[str]]
            - options : Optional[Dict[str, Any]]
    script : str
        The script source code.
    language : ScriptLanguage
        The programming language used for the script.

    Examples
    --------
    Creating and serializing a Python script function:

    >>> script_code = "print('Hello from Datagrok!')"
    >>> s = Script(script=script_code, language=ScriptLanguage.Python, name="Hello Script")
    >>> s.to_dict()
    ... {
    ...     "id": None,
    ...     "name": "Hello Script",
    ...     "friendlyName": None,
    ...     "createdOn": None,
    ...     "updatedOn": None,
    ...     "source": "script",
    ...    "description": None,
    ...     "params": [],
    ...     "namespace": None,
    ...     "tags": [],
    ...     "options": {},
    ...     "script": "print('Hello from Datagrok!')",
    ...     "language": "python"
    ... }
    """
    SOURCE = "script"

    def __init__(self, script: str, language: Optional[ScriptLanguage]=ScriptLanguage.Python, **kwargs):
        kwargs.pop("source", None)
        super().__init__(source=Script.SOURCE, **kwargs)
        self.script = script
        self.language = language

    def to_dict(self):
        d = super().to_dict()
        d["script"] = self.script
        d["language"] = self.language.value
        return d
    
    @classmethod
    def _from_dict(cls, data: dict) -> "Script":
        return cls(
            script=data.get("script"),
            language=ScriptLanguage(data.get("language")),
            id=data.get("id"),
            name=data.get("name"),
            friendly_name=data.get("friendlyName"),
            description=data.get("description"),
            created_on=datetime.fromisoformat(data["createdOn"]) if data.get("createdOn") else None,
            updated_on=datetime.fromisoformat(data["updatedOn"]) if data.get("updatedOn") else None,
            source=data.get("source"),
            params=[FuncParam.from_dict(d) for d in data.get("params")] if data.get("params") else [],
            namespace=data.get("namespace"),
            tags=data.get('tags', []),
            options=data.get("options")
        )
    