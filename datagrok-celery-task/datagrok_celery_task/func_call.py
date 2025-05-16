from enum import Enum
from typing import Union, Any, Dict, List


class Type(str, Enum):
    INT = 'int'
    STRING = 'string'
    BOOL = 'bool'
    FLOAT = 'double'
    BIG_INT = 'bigint'
    DATA_FRAME = 'dataframe'
    DATE_TIME = 'datetime'
    FILE = 'file'
    BLOB = 'blob'


class FuncCallStatus(str, Enum):
    RUNNING = 'Running'
    QUEUED = 'Queued'
    COMPLETED = 'Completed'
    ERROR = 'Error'
    CANCELED = 'Canceled'
    CANCELING = 'Canceling'
   

class FuncCallParam:
    def __init__(self, **kwargs):
        self.name = kwargs.get("name")
        self.property_type = Type(kwargs.get("property_type"))
        self.value = kwargs.get("value")
        self.is_input = kwargs.get("is_input")

    @property
    def is_streamable(self) -> bool:
        return self.property_type == Type.DATA_FRAME or self.property_type == Type.BLOB or self.property_type == Type.FILE


class FuncCall:
    default_binary_batch_size: int = 2048000

    def __init__(self, call: Dict[str, Any]):
        self.id = call.get("id")
        self.func: Dict[str, Any] = call.get("func", {})
        self.options: Dict[str, Any] = call.get("options", {})
        self.aux: Dict[str, Any] = call.get("aux", {})
        self.status = FuncCallStatus.RUNNING

        parameter_values: Dict[str, Any] = call.get("parameterValues", {})
        func_params: List[Dict] = self.func.get("params") or []
        self.params: List[FuncCallParam] = []

        for param in func_params:
            param_name = param.get("name")
            value = parameter_values.get(param_name, None)
            self.params.append(FuncCallParam(name=param_name, value=value, property_type=param.get("propertyType"), is_input=param.get("isInput", True)))

        self.requires_pipe = any(param.is_streamable for param in self.params)
        self.error_message = None
        self.error_stack_trace = None

    @property
    def use_parquet_transfer(self) -> bool:
        return self.options.get("isParquet", False)

    @property
    def user_api_key(self) -> Union[str, None]:
        return self.aux.get("USER_API_KEY", None)
    
    @property
    def binary_batch_size(self) -> int:
        return self.aux.get("batchSize", self.default_binary_batch_size)

    @property
    def input_params(self):
        return [x for x in self.params if x.is_input]  

    @property
    def output_params(self):
        return [x for x in self.params if not x.is_input]  

    def to_json(self):
        return {
            "id": self.id,
            "parameterValues": {p.name: p.value for p in self.params},
            "func": self.func,
            "options": self.options,
            "aux": self.aux,
            "status": self.status.value,
            "errorMessage": self.error_message,
            "errorStackTrace": self.error_stack_trace
        }
    