import re
from enum import Enum
import io
import pandas as pd
import numpy as np
from datetime import datetime
from typing import Any

import pyarrow.parquet as pq
import pyarrow as pa
import json

from .func_call import Type, FuncCallParam


class DatagrokFanoutType(str, Enum):
    CALL = "call"
    LOG = "log"
    ACCEPTED = "accepted"
    PROGRESS = "progress"


def camel_to_snake(name: str) -> str:
    s1 = re.sub(r'(.)([A-Z][a-z]+)', r'\1_\2', name)
    s2 = re.sub(r'([a-z0-9])([A-Z])', r'\1_\2', s1)
    return s2.lower()


class InputValueProcessor:
    @staticmethod
    def set_datetime(param: FuncCallParam):
        if not isinstance(param.value, (str, datetime)):
            raise ValueError(f"Incorrect input type for {param.name}. Expected str or datetime, got {type(param.value)}")
        
        param.value = datetime.strptime(param.value, "%Y-%m-%dT%H:%M:%S.%f%z") if type(param.value) == str else param.value

    @staticmethod
    def set_dataframe(param: FuncCallParam, use_parquet: bool):
        if not isinstance(param.value, (bytes, bytearray)):
            raise ValueError(f"Incorrect input type for {param.name}. Expected bytes or bytearray, got {type(param.value)}")
        if use_parquet:
            buffer = pa.BufferReader(param.value)
            param.value = pq.read_table(buffer).to_pandas()
        else:    
            param.value = pd.read_csv(io.StringIO(param.value.decode("utf-8")), skip_blank_lines=False)

    @staticmethod
    def set_bigint(param: FuncCallParam):
        if not isinstance(param.value, (str, int)):
            raise ValueError(f"Incorrect input type for {param.name}. Expected str or int, got {type(param.value)}")
        try:
            param.value = int(param.value)
        except (ValueError, TypeError):
            raise ValueError(f"Incorrect return type for {param.name}. Expected str-representable bigint or int.")

    type_map = {
        Type.BIG_INT: set_bigint.__func__,
        Type.DATE_TIME: set_datetime.__func__,
        Type.DATA_FRAME: set_dataframe.__func__,
    }


class ReturnValueProcessor:
    @staticmethod
    def set_int(param: FuncCallParam, val: Any):
        try:
            param.value = int(val)
        except (ValueError, TypeError):
            raise ValueError(f"Incorrect return type for {param.name}. Expected int, got {type(val)}")

    @staticmethod
    def set_float(param: FuncCallParam, val: Any):
        try:
            param.value = float(val)
        except (ValueError, TypeError):
            raise ValueError(f"Incorrect return type for {param.name}. Expected float, got {type(val)}")

    @staticmethod
    def set_bool(param: FuncCallParam, val: Any):
        try:
            param.value = bool(val)
        except (ValueError, TypeError):
            raise ValueError(f"Incorrect return type for {param.name}. Expected bool, got {type(val)}")

    @staticmethod
    def set_str(param: FuncCallParam, val: Any):
        try:
            param.value = str(val)
        except (ValueError, TypeError):
            raise ValueError(f"Incorrect return type for {param.name}. Expected str, got {type(val)}")

    @staticmethod
    def set_bigint(param: FuncCallParam, val: Any):
        try:
            param.value = str(val)
        except (ValueError, TypeError):
            raise ValueError(f"Incorrect return type for {param.name}. Expected str-representable bigint.")

    @staticmethod
    def set_datetime(param: FuncCallParam, val: Any):
        if not isinstance(val, datetime):
            raise ValueError(f"Incorrect return type for {param.name}. Expected datetime, got {type(val)}")
        param.value = val.isoformat()

    @staticmethod
    def set_blob(param: FuncCallParam, val: Any):
        if not isinstance(val, (bytes, bytearray)):
            raise ValueError(f"Incorrect return type for {param.name}. Expected bytes, got {type(val)}")
        param.value = val

    @staticmethod
    def set_dataframe(param: FuncCallParam, val: Any, use_parquet: bool):
        if not isinstance(val, pd.DataFrame):
            raise ValueError(f"Incorrect return type for {param.name}. Expected pandas.DataFrame, got {type(val)}")
        if use_parquet:
            ReturnValueProcessor.flatten_columns(val)
            ReturnValueProcessor.fill_nulls_for_export(val)
            ReturnValueProcessor.downcast_int64_to_int32(val)
            buffer = io.BytesIO()
            val.to_parquet(buffer, engine="pyarrow", compression="snappy")
            buffer.seek(0)
            param.value = buffer.getvalue()
        else:
            param.value = val.to_csv(index=False).encode("utf-8")

    @staticmethod
    def flatten_columns(df, sep="."):
        c = df.columns
        if isinstance(c, pd.MultiIndex):
            df.columns = [sep.join(map(str, filter(None, x))) for x in c.values]
        elif isinstance(c, pd.CategoricalIndex):
            df.columns = c.astype(str)
        elif not all(isinstance(x, str) for x in c):
            df.columns = list(map(str, c))

        for col in df.columns:
            s = df[col].dropna().head(10).astype(object)
            if not s.empty and s.apply(lambda x: isinstance(x, (list, dict))).any():
                df[col] = df[col].apply(lambda x: json.dumps(x) if pd.notnull(x) else "")

    @staticmethod
    def fill_nulls_for_export(df):
        for c in df.columns:
            s = df[c]
            if pd.api.types.is_integer_dtype(s):
                df[c] = s.fillna(-2147483648)
            elif pd.api.types.is_float_dtype(s):
                df[c] = s.fillna(2.6789344063684636e-34)
            elif pd.api.types.is_object_dtype(s) or pd.api.types.is_categorical_dtype(s):
                df[c] = s.fillna("")

    @staticmethod
    def downcast_int64_to_int32(df):
        def is_nullable_integer_dtype(x):
            return (
                pd.api.types.is_extension_array_dtype(x)
                and pd.api.types.pandas_dtype(x).name.startswith("Int")
            )

        mn, mx = -2_147_483_648, 2_147_483_647
        for c in df.columns:
            s = df[c]
            if pd.api.types.is_integer_dtype(s):
                v = s.dropna().astype(np.int64)
                if not v.empty and v.min() >= mn and v.max() <= mx:
                    df[c] = (
                        s.astype("Int32")
                        if is_nullable_integer_dtype(s)
                        else s.astype("int32")
                    )

    type_map = {
        Type.INT: set_int.__func__,
        Type.FLOAT: set_float.__func__,
        Type.BOOL: set_bool.__func__,
        Type.STRING: set_str.__func__,
        Type.BIG_INT: set_bigint.__func__,
        Type.DATE_TIME: set_datetime.__func__,
        Type.BLOB: set_blob.__func__,
        Type.DATA_FRAME: set_dataframe.__func__,
    }
