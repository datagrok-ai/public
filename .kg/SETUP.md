# Knowledge graph — setup

The graph lives at `.kg/` and ships with the public submodule. To
run `qq.py`, `build.py`, the enrichers, the visualization, or the agent
harness in `.claude/`, you need a Python venv with three packages
installed.

## TL;DR — one-shot bootstrap

```powershell
# Windows PowerShell (run from the public root)
.kg\scripts\bootstrap.ps1
```

```bash
# macOS / Linux / Git Bash
bash .kg/scripts/bootstrap.sh
```

Both scripts are idempotent. They:
1. Create `.kg/.venv/` if missing
2. Install `requirements.txt` into it
3. Run a smoke test: `qq.py "MATCH (n) RETURN count(n) LIMIT 1"`
4. Print the venv Python path so you can copy it into your `PATH` or
   alias if you want.

If the smoke test fails because the Kuzu DB is missing (fresh clone),
the script will offer to run `build.py` to materialize it from JSONL.
That takes ~3 minutes the first time.

## Manual setup (if you don't trust the script)

```powershell
# 1. Create venv
py -m venv .kg/.venv

# 2. Install requirements
.kg/.venv/Scripts/python.exe -m pip install -r .kg/requirements.txt

# 3. Materialize Kuzu DB from canonical JSONL (only needed on fresh clone)
.kg/.venv/Scripts/python.exe .kg/build.py

# 4. Sanity check
.kg/.venv/Scripts/python.exe .kg/qq.py "MATCH (n) RETURN count(n)"
```

`bash` users: replace `.kg/.venv/Scripts/python.exe` with
`.kg/.venv/bin/python` and use `python3 -m venv` if `py` isn't
on your PATH.

## Daily use

After bootstrap, you have two convenient front-doors:

```powershell
# Persistent server (recommended for any iterative work)
.kg/.venv/Scripts/python.exe .kg/qq.py "<your cypher>"

# One-shot CLI (no server)
.kg/.venv/Scripts/python.exe .kg/query.py --demos
```

The `qq.py` server starts automatically on first call and runs in the
background; subsequent calls hit it in ~5–50 ms. Stop it with
`qq.py --stop` (do this before any rebuild — the server holds the DB
lock).

## What gets installed

`requirements.txt` is intentionally minimal:

| Package | Why |
|---|---|
| `pydantic >= 2.6` | Schema models for entities/edges |
| `kuzu >= 0.6` | Embedded graph DB (MIT-licensed, single Python wheel) |
| `pyyaml >= 6.0` | Reading `environments/*.yaml` for ScriptEnvironment extraction |
| `tomli >= 2.0` (Py < 3.11 only) | Reading some package configs |

No system dependencies. No Docker. No cloud. The whole graph lives in
`.kg/data/*.jsonl` (~12 MB, committed) and `.kg/kg.kuzu/`
(~400 MB, gitignored, rebuildable).

## What if something goes wrong

| Symptom | Fix |
|---|---|
| `ModuleNotFoundError: kuzu` | You're using system `py` instead of the venv. Use the absolute venv path. |
| `IO exception: Could not set lock on file` | The `qq.py` server is running. `qq.py --stop` then retry. |
| Cold start > 5s | First `qq.py` call spawns the server — that's the 0.8s warmup, then queries are 5–50ms. |
| `MemoryError` during build | A previous extractor run left huge JSONL artifacts. `build.py --clean` then re-apply enrichers (see `SESSION_NOTES.md`). |
| Empty result on a query you expected to work | The graph is a snapshot. After local code changes, run `build.py --packages <YourPkg>` to refresh just that package. |

For deeper context (architecture, lessons learned, open issues), read
[`SESSION_NOTES.md`](SESSION_NOTES.md).
