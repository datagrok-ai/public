# Scripting

Per-language script execution containers. This package ships **no client code** — it exists
so the platform manages the five kernel gateways as on-demand Docker containers instead of
always-on infrastructure services.

| Container entity | Image | CPU / RAM |
|---|---|---|
| `scripting-python` | `datagrok/jkg_python` | 2 / 8 GB |
| `scripting-r` | `datagrok/jkg_r` | 1 / 4 GB |
| `scripting-octave` | `datagrok/jkg_octave` | 1 / 2 GB |
| `scripting-julia` | `datagrok/jkg_julia` | 1 / 4 GB |
| `scripting-nodejs` | `datagrok/jkg_nodejs` | 1 / 2 GB |

Entity names follow the platform rule `kebab(package) + '-' + <dockerfiles folder>`, and the
server derives them from the script's language — see `ScriptingService.ensureLanguageWorker`.

## How a script reaches its container

1. A script runs; `Func.queueInterceptor` sees it needs a capability this node lacks.
2. Before the call is queued, the interceptor resolves `scripting-<language>`, starts it if
   it is stopped (`on_demand`), and bumps `last_active`.
3. It then blocks until that container's worker is actually **consuming**
   `scripting.jupyter.<language>` — a started container is not the same as a subscribed one,
   and publishing early would drop the message.
4. Only then is the call published.

The wait budget is `queueSettings.containerWakeTimeout` (default 10 min), deliberately
independent of `taskPickupTimeout` (60 s): a cold `jkg_python` has to pull ~7.7 GB and boot
conda before it can subscribe.

## Setup wizard

The initial setup wizard has a **Scripting** page (right after the health check). JavaScript and
Grok are shown checked and disabled — they need no container — and the five container languages
are offered with their download sizes, Python pre-checked. Picking any of those posts to
`/scripts/languages/setup`, which installs this package and starts the chosen containers.

That call returns as soon as the work is scheduled — a cold `jkg_python` pull outlasts any
sensible request timeout — and reports each step back over `Events.CLIENT_PUSH` as a
`ServerTaskProgress`. The wizard renders those pushes inline and balloons them, so the user is
told *why* the wait is long instead of watching a spinner. Implementation:
`ScriptingService.setupLanguages`.

## Idle shutdown

`shutdown_timeout: 2880` (48 h) with `on_demand: true`. Datlas's reaper stops a container
whose `last_active` is older than that; `on_demand` is what allows it to be restarted. The
two must always be set together — `shutdown_timeout` alone stops a container that nothing
will ever bring back.

`last_active` is normally only bumped by HTTP proxy traffic, which these containers never
receive, so the interceptor bumps it on every dispatch. Without that the 48 h would run from
container start rather than from last use.

## Images

There are no Dockerfiles here — each `dockerfiles/<lang>/container.json` names the published
image to run, and `grok publish` builds and pushes nothing. The images are built by CI from
`deploy/jkg_<lang>/`, because they need the core source tree to compile the datlas worker
bundle, which a package build context does not have. To roll a language forward, bump the
tag in its `container.json`.

The first time the spawner resolves one of these tags to Docker Hub it also copies it into the
cluster's internal registry, in the background (`spawn/registry_mirror.py`). Every later pull —
including the one after the 48 h idle shutdown, and the ones on other nodes — is then a
cluster-local transfer instead of another multi-GB trip to the Hub. Sizes below are the layer
sums, which is what actually transfers.

Design, measurements and the vulnerability scan: `core/docs/features/jkg-language-split/`.
