# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is the **Datagrok Python API** (`datagrok-api`) — a Python client library for the Datagrok platform's REST API. It wraps the public API endpoints (`/public/v1/*`) and uses pandas DataFrames for table data exchange.

## Commands

```bash
# Install (editable, with test deps)
uv pip install -e ".[test]"

# Run all tests (requires a running Datagrok server)
DATAGROK_API_URL=http://localhost:8082 DATAGROK_API_KEY=unit_test_token uv run pytest

# Or pass via CLI options
uv run pytest --url http://localhost:8082 --key unit_test_token

# Run a single test file
uv run pytest tests/test_tables_client.py --url http://localhost:8082 --key unit_test_token

# Run a specific test
uv run pytest tests/test_tables_client.py::test_upload_and_download --url http://localhost:8082 --key unit_test_token

# Build distribution
uv build

# Generate docs (outputs Markdown for Docusaurus)
uv run python generate_docs.py
```

No linting configuration exists in this repo.

## Architecture

Three-layer design, all wired together by `DatagrokClient` (facade):

```
DatagrokClient          ← Entry point, context manager, aggregates all resource clients
├── resources/          ← Resource clients (one per API domain)
│   ├── tables.py       ← TablesClient: upload/download DataFrames (CSV over HTTP)
│   ├── files.py        ← FilesClient: upload/download/sync files, uses .sync_meta.json for eTags
│   ├── functions.py    ← FunctionsClient: call functions, create scripts/queries
│   ├── users.py        ← UsersClient: current user, find users
│   ├── groups.py       ← GroupsClient: group hierarchy, membership management
│   ├── connections.py  ← ConnectionsClient: CRUD for data connections, test connectivity
│   └── shares.py       ← SharesClient: share entities with groups/users
├── models/             ← Data models with to_dict()/from_dict() serialization
│   ├── model.py        ← Base classes: Model (has id), NamedModel (adds namespace, name, timestamps)
│   ├── user.py         ← User (status enum, login, group membership)
│   ├── group.py        ← Group + GroupRelation (parent-child with admin flags)
│   ├── data_connection.py  ← DataConnection + Credentials (dynamic attrs via __getattr__)
│   ├── func.py         ← Func base + FuncParam + PropertyType enum
│   ├── script.py       ← Script subclass (registered in Func._registry)
│   ├── query.py        ← DataQuery subclass (registered in Func._registry)
│   └── share_response.py
└── http_client.py      ← Thin wrapper around requests.Session with auth header
```

## Key Patterns

- **Polymorphic deserialization**: `Func.from_dict()` uses a `_registry` dict keyed by `source` field to instantiate the correct subclass (`Script` for "script", `DataQuery` for "data-query").
- **Dynamic attributes on connections**: `DataConnection` and `Credentials` use `__getattr__`/`__setattr__` to store arbitrary key-value pairs in `self.parameters`, allowing `conn.bucket` or `creds.accessKey` without explicit fields.
- **Grok name format**: Entity names use colon-separated format (`"Namespace:Name"`) which gets converted to dot notation for URL paths (`"Namespace.Name"`).
- **Tests run against a real server**: All tests require a running Datagrok instance. Configure via `DATAGROK_API_URL`/`DATAGROK_API_KEY` env vars or `--url`/`--key` pytest options. Tests that create entities clean up after themselves.
