# Agentic force for scientists

Grokky turns Claude into a context-aware agent that knows what the platform can do, what the
company's procedures are, and what the user has taught it. This is powered by the **skills &
knowledge system** — a sync layer that keeps the Claude runtime's local filesystem in sync with
files from three sources: the user, published packages, and shared connections.

---

## Skills overview

A "skill" is any file (typically `.md`) that gets included in Claude's context at session start.
Skills teach Claude about available functions, domain procedures, infrastructure, or anything else
that makes it more effective. Three scopes are supported:

| Scope | Source | Example |
|-------|--------|---------|
| **Personal** | User's `My Files/agents/` folder | "Always use project X defaults when running ADME" |
| **Package** | `agents/` folder inside a published package | "For chemical space, call the ChemicalSpace function" |
| **Shared** | Shared `agents/` folders from other users' My Files | Company SOPs, team-specific procedures |

**Package skills** are the primary mechanism for plugins to describe their capabilities.
Standard platform packages (like Chem) ship with an `agents/` folder — these live in the
public workspace alongside the package source. User-created packages also support `agents/`
folders, but since they are personal, their skills sync into the user's own
`/users/{userId}/agents/` directory. A Chem plugin might ship:

```
for chemical space, call the ChemicalSpace function
```

**Shared skills** enable company-specific procedures or instructions:

```
When registering a new chemical compound, do the following
- make sure it's a valid structure
- ask user to associate it with a project
- register in Global Compound Registration System
- kick off ADME prediction
- email John Smith the link to ADME results
```

This becomes both a mechanism that defines what is available to the user (via plugin privileges)
and a mechanism for efficient working.

---

## Folder creation & sharing

On init, Grokky creates `agents/` as a directory inside My Files (with a `README.md` placeholder).
Users place their personal knowledge files here.

To share skills with others, a user shares the `agents/` folder through the platform UI (right-click
→ Share). This creates a sub-connection that other users can access. The sync
layer discovers these shared connections by querying `/public/v1/connections?text=agents` and
filtering for connections named `MyFilesAgents`. No tags or special setup required — sharing the
folder is the only step.

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

---

## Sync

All three scopes sync independently using timestamp-based incremental sync — only new or changed
files are downloaded, and files deleted on the remote side are cleaned up locally. Naming
conventions prevent collisions across sources:

| Scope | Local naming | Trigger |
|-------|-------------|---------|
| Personal | `agents/{filename}` | File events (`d4-file-event`, `onFileEdited`) |
| Package | `agents/{PackageName}-{filename}` | `onPackageLoaded` event + 15-min poll |
| Shared | `agents/shared-{connId}-{filename}` | 15-minute polling interval |

Package sync uses timestamp caching (`updatedOn`) to skip unchanged packages — ZIP downloads
only happen when a package has actually been updated.

On the first user message, all three scopes sync (`scope=all`). After that, event-driven syncs
use a narrow scope so only the relevant category is re-checked. A concurrency guard prevents
duplicate syncs for the same user.

---

## Open design questions

1. **File count scaling** — the system prompt currently lists up to 50 agent files. As the
   number of packages and shared connections grows, a smarter selection or indexing strategy
   may be needed (e.g., relevance-based filtering, embedding search).

2. **Conflict resolution** — if a personal file and a package file have the same effective name,
   the current naming convention prevents filesystem collisions, but Claude sees both. A priority
   model (personal > shared > package) could help resolve contradictions in instructions.