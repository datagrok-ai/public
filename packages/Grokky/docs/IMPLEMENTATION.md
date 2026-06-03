# Skills & knowledge sync — implementation

Technical implementation of the skills & knowledge sync system. See [VISION.md](VISION.md)
for the user-facing concept.

---

## Folder bootstrap

On init, Grokky creates an `agents/` directory in the user's My Files with a `README.md`
placeholder.

---

## Shared folder discovery

When a user shares their `agents/` folder, the platform creates a sub-connection named
`MyFilesAgents`. The sync layer discovers these by querying
`/public/v1/connections?text=agents` and filtering by name. No tags or registration are
required — the named sub-connection is the entire contract.

---

## Working directory

The Claude runtime maintains a per-user directory on the container filesystem:

```
/users/{userId}/
├── agents/                          ← all knowledge files land here
│   ├── my-notes.md                  ← personal (from My Files/agents/)
│   ├── <UserPackageName>-skills.md  ← package (from <UserPackageName> package's agents/)
│   └── shared-abc123-onboard.md     ← shared (from another user's agents/)
└── workspace → /workspace           ← symlink to the cloned Datagrok public repo
```

The user ID is extracted from the JWT payload (`sub` or `usr.id`). The `workspace` symlink gives
Claude read access to the full JS API source, packages, samples, and documentation.

The starting point for Claude is determined by whether any user knowledge exists — if the user has
agent files, Claude's working directory is the user directory; otherwise it defaults to `/workspace`.

All agent files are listed in Claude's system prompt (up to 50) so it knows what knowledge is
available without needing to scan the filesystem.

User-created packages sync into `/users/{userId}/agents/` rather than the public workspace,
since their skills are personal.

---

## Sync

All three scopes sync independently using incremental sync — only new or changed files are
downloaded, and files deleted on the remote side are cleaned up locally. Naming conventions
prevent collisions across sources:

| Scope | Source | Local naming | Trigger |
|-------|--------|--------------|---------|
| Personal | User's `My Files/agents/` folder | `agents/{filename}` | File events (`d4-file-event`, `onFileEdited`) |
| Package | `agents/` folder inside a published package | `agents/{PackageName}-{filename}` | `onPackageLoaded` event + on-demand TTL refresh |
| Shared | `MyFilesAgents` sub-connections from other users | `agents/shared-{connId}-{filename}` | On-demand TTL refresh |

Package sync keys on the package's `buildHash` to skip unchanged packages. When a build changes,
the runtime lists the published-files folder
(`/packages/published/files/{name}/{version}/{buildHash}/{buildNumber}/agents/`) and downloads
each agent file individually. The per-user cache lives in process memory only and tracks
`buildHash` plus the files written for that build, so the previous build's files can be removed
before the new ones land.

Shared sync uses `updatedOn` timestamp caching. The seen version is recorded whether or not the
download succeeds, so a file that can't be fetched (e.g. invalid JSON) is not retried on every
refresh — it retries only when its `updatedOn` changes.

The first user message syncs all scopes blockingly (`scope=all`); later `scope=all` turns serve
on-disk state immediately and refresh the package and shared scopes in the background once a
15-minute TTL elapses, so a turn never pays sync latency. A concurrency guard (`syncInFlight` per
user) prevents duplicate syncs.

---

## Package knowledge index

Packages can ship a structured `agents/package-knowledge.yaml` file that describes the package's
capabilities in a machine-readable format:

```yaml
packageName: Admetica
description: ADMET property prediction for small molecules using ChemProp v2 neural networks
keywords:
  - admet
  - absorption
  - molecules
  - smiles
  - chemprop
overview: |
  Predicts ADMET molecular properties using ChemProp v2 neural networks.

  Single molecule: click on a molecule cell — the Admetica panel appears
  in the context panel (Biology section) with predictions.
apiRef: src/package-api.ts
docsRef: ../../help/datagrok/solutions/domains/chem/chem.md
```

At startup, the runtime scans `/workspace/packages/*/agents/package-knowledge.yaml` on the local
filesystem (no network calls, no ZIP downloads) and aggregates all found files into a unified
markdown table. This table is injected directly into Claude's system prompt as
`## Available Packages`, giving Claude immediate awareness of all package capabilities without
requiring a tool call.

The index is cached after the first scan — subsequent messages reuse the cached result.

Claude is instructed to consult this table first when a user asks about platform capabilities,
and to read the full `package-knowledge.yaml` for details before falling back to code search.

**Future direction:** consolidate the package index, agent file listings, and other dynamic context
into a generated `CLAUDE.md` or skills file rather than appending sections to the system prompt.
This would align with how Claude Code natively consumes project context and make the knowledge
layer easier to inspect and debug.

---

## Open design questions

1. **File count scaling** — the system prompt currently lists up to 50 agent files. As the
   number of packages and shared connections grows, a smarter selection or indexing strategy
   may be needed (e.g., relevance-based filtering, embedding search).

2. **Conflict resolution** — if a personal file and a package file have the same effective name,
   the current naming convention prevents filesystem collisions, but Claude sees both. A priority
   model (personal > shared > package) could help resolve contradictions in instructions.
