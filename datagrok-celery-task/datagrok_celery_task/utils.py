import re
from enum import Enum
import io
import pandas as pd
from datetime import datetime
from typing import Any

import pyarrow.parquet as pq
import pyarrow as pa

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
        cols = df.columns

        if isinstance(cols, pd.MultiIndex):
            df.columns = [sep.join(map(str, filter(None, x))) for x in cols.values]
        elif isinstance(cols, pd.CategoricalIndex):
            df.columns = cols.astype(str)
        elif not all(isinstance(x, str) for x in cols):
            df.columns = list(map(str, cols))

        for col in df.columns:
            sample = df[col].dropna().head(10)
            if not sample.empty and sample.apply(lambda x: isinstance(x, (list, dict))).any():
                df[col] = df[col].apply(lambda x: json.dumps(x) if pd.notnull(x) else "")

    @staticmethod
    def fill_nulls_for_export(df):
        for col in df.columns:
            series = df[col]
            if pd.api.types.is_integer_dtype(series):
                df[col] = series.fillna(-2_147_483_648)
            elif pd.api.types.is_float_dtype(series):
                df[col] = series.fillna(2.6789344063684636e-34)
            elif (
                pd.api.types.is_object_dtype(series)
                or pd.api.types.is_string_dtype(series)
                or pd.api.types.is_categorical_dtype(series)
            ):
                df[col] = series.fillna("")

    @staticmethod
    def downcast_int64_to_int32(df):
        INT32_MIN, INT32_MAX = -2_147_483_648, 2_147_483_647

        for col in df.columns:
            series = df[col]
            if pd.api.types.is_integer_dtype(series):
                values = series.dropna().astype(np.int64)
                if not values.empty and values.min() >= INT32_MIN and values.max() <= INT32_MAX:
                    df[col] = (
                        series.astype("Int32")
                        if pd.api.types.is_nullable_integer_dtype(series)
                        else series.astype("int32")
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
