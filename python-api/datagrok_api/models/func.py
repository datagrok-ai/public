from datetime import datetime
from enum import Enum
from typing import Optional, Any, Dict, List, Type
from datagrok_api.models.model import NamedModel

class PropertyType(str, Enum):
    INT = 'int'
    STRING = 'string'
    BOOL = 'bool'
    FLOAT = 'double'
    BIG_INT = 'bigint'
    DATA_FRAME = 'dataframe'
    DATE_TIME = 'datetime'
    FILE = 'file'
    BLOB = 'blob'
    COLUMN = 'column'
    COLUMN_LIST = 'column_list'
    GRAPHICS = 'graphics'
    DYNAMIC = 'dynamic'

class FuncParam:
    def __init__(self, name: str, type: PropertyType, is_input: bool = True):
        self.name = name
        self.type = type
        self.is_input = is_input

    def to_dict(self):
        return {
            "name": self.name,
            "propertyType": self.type,
            "isInput": self.is_input
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> 'FuncParam':
        return cls(name=data.get("name"), type=PropertyType(data.get("propertyType")), is_input=data.get("isInput", True))

class Func(NamedModel):
    """
    Represents a function in the Datagrok platform.

    In Grok, nearly every executable action—whether it's a data query, 
    a script, a UI action, or a predictive model—is represented as a 
    `Func`. This unified abstraction provides a consistent way to 
    define, discover, secure, and execute platform capabilities.

    Examples of functions in Grok include:
        * Querying an external PostgreSQL database
        * Executing JavaScript code in the browser that uses Grok APIs
        * Performing mathematical calculations (e.g., `Sin(PI)`)
        * Modifying datasets (e.g., deleting a column)
        * Sending an email
        * Applying a predictive model to a dataset
        * Calculating molecular properties via a Python script
        * Displaying interactive dialogs

    Despite their differences (server-side vs client-side execution, 
    computational vs UI-based tasks), all functions share a common 
    set of capabilities:
        * **Scriptable**: Callable from the console or scripts.
        * **Findable**: Searchable through the "Help | Functions" interface.
        * **Introspectable**: Metadata about parameters is programmatically accessible.
        * **Secure**: Access can be restricted via privileges or user groups.
        * **Auditable**: Execution history and parameters are tracked.
        * **Runnable**: The platform can dynamically generate UI for parameters.
        * **Linkable**: Functions can be linked in dashboards, conversations, etc.
        * **Composable**: Usable in workflow designers and query transformations.

    Attributes
    ----------
    source : Optional[str], default='function'
        Identifies the origin or type of function (e.g., `"function"`, `"sql-query"`, `"js-script"`).
    options : Optional[Dict[str, Any]], default=None
        Arbitrary function-specific configuration options.
    tags : Optional[List[str]], default=None
        Keywords or categories associated with the function.
    description : Optional[str], default=None
        A human-readable description of the function's purpose.
    params : Optional[List[FuncParam]], default=None
        List of parameter definitions for the function.
    **kwargs
        Additional arguments passed to the parent `NamedModel` initializer,
        such as `id`, `name`, `friendly_name`, `created_on`, `updated_on`, and `namespace`.
    params : list of FuncParam
        Parameter definitions (both input and output).
    source : str
        Function type or origin.
    tags : list of str
        Function categories or keywords.
    description : str
        Human-readable description.
    options : dict
        Function-specific settings.
    input_params : list of FuncParam
        Parameters marked as inputs to the function.
    output_params : list of FuncParam
        Parameters marked as outputs from the function.
    """
    SOURCE = "function"
    _registry: Dict[str, Type["Func"]] = {}

    def __init__(self, source: Optional[str] = 'function', options: Optional[Dict[str, Any]] = None, tags: Optional[List[str]] = None, description: Optional[str]=None, params: Optional[List[FuncParam]] = None, **kwargs):
        super().__init__(**kwargs)
        self.params = params or []
        self.source = source
        self.tags = tags or []
        self.description = description
        self.options = options or dict()

    @property
    def input_params(self):
        return [x for x in self.params if x.is_input]  

    @property
    def output_params(self):
        return [x for x in self.params if not x.is_input]      

    def to_dict(self):
        return {
            "id": self.id,
            "name": self.name,
            'friendlyName': self.friendly_name,
            "createdOn": self.created_on.isoformat() if self.created_on else None,
            "updatedOn": self.updated_on.isoformat() if self.updated_on else None,
            "source": self.source,
            "description": self.description,
            "params": [param.to_dict() for param in self.params],
            "namespace": self.namespace,
            "tags": self.tags,
            "options": self.options
        }
    
    @classmethod
    def register_subclass(cls, source: str):
        def decorator(subclass):
            Func._registry[source] = subclass
            return subclass
        return decorator   
    
    @classmethod
    def from_dict(cls, data: dict) -> "Func":
        source = data.get("source")
        subclass = Func._registry.get(source, cls)
        return subclass._from_dict(data)
    
    @classmethod
    def _from_dict(cls, data: dict) -> "Func":
        return cls(
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
    