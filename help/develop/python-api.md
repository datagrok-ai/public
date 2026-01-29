---
title: Python API
---

The datagrok-api Python library is designed for integration with the Datagrok platform via REST API. It provides a convenient 
interface for working with data, tables, files, and other Datagrok objects.

## Installation

The library can be installed using pip:

```shell
pip install datagrok-api
```

## Usage

```python
from datagrok_api import DatagrokClient

# Initialize with your server URL and API token
with DatagrokClient("https://public.datagrok.ai/api", "Bearer <your-token>") as grok:
    me = grok.users.current()
    print(f"Logged in as: {me.first_name} {me.last_name}")
```

## API structure

The Python API is built around a central `DatagrokClient` and a set of dedicated resource clients. `DatagrokClient` is the main entry point to the platform: it handles authentication 
and HTTP communication and exposes resource-specific clients as its properties.

Each resource module implements a dedicated client responsible for a specific domain, including user-related operations (`users`), 
group search and management (`groups`), upload and download of tables as Pandas DataFrames (`tables`), file upload and download (`files`), 
execution of server-side functions (`functions`), data connection management (`connections`), and shared entities and permissions (`shares`).

Typical usage follows this pattern:

```python
grok.users.current()
grok.tables.download("JohnDoe:MyTable")
grok.groups.find("Chemists")
```

The models package contains typed representations of Datagrok entities (such as User, Group, DataConnection, Query, and Script) that are used as inputs 
and return values by the resource clients. All models share a common base class defined in model.py.
