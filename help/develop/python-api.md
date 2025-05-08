# Python API

The datagrok-api Python library is designed for integration with the Datagrok platform via REST API. It provides a convenient 
interface for working with data, tables, files, and other Datagrok objects.

## Installation

The library can be installed using pip:

```bash
pip install datagrok-api
```

The main class of the library is `DatagrokClient`, which provides methods for interacting with the Datagrok API.

```python
from datagrok_api import DatagrokClient
api = DatagrokClient('your_api_key', 'datagrok_api_url')****
```

## API structure

API contains several modules:

* [**datagrok_client**](https://datagrok.ai/api/py/datagrok_client) is the main module that allows interaction with Datagrok platform,
* [**group**](https://datagrok.ai/api/py/group) contains classes for groups management,
* [**model**](https://datagrok.ai/api/py/model) contains single class, base class for [**group**](https://datagrok.ai/api/py/group) classes