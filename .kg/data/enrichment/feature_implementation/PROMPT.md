# Feature implementation/test mapping

You are mapping each Feature in a single Datagrok plugin to the source
files that **implement** it and the source files that **test** it.

## Why this exists

Datagrok plugins follow a wrapper pattern:
- `src/package.g.ts` is auto-generated and just contains `//name:` annotation
  blocks + thin `export async function foo() { return PackageFunctions.foo() }`
  wrappers.
- `src/package.ts` declares `class PackageFunctions` whose decorated `static`
  methods call into the real implementation in `src/widgets/`, `src/utils/`,
  `src/viewers/`, `src/analysis/` etc.

A naive walker only sees `package.g.ts`. We have a deterministic AST walker
that follows the `package.ts` import chain — that result is in each feature's
`ast_baseline.implements`. Your job is to **extend and correct** that baseline.

## Your input

A JSON object:

```
{
  "package": { "id", "name", "friendly_name", "category", "description" },
  "files": {
    "src":     [paths under <pkg>/src/ excluding tests/, .g.ts, package-api.ts],
    "tests":   [paths under <pkg>/src/tests/ + package-test.ts],
    "scripts": [paths under <pkg>/scripts/],
    "queries": [paths under <pkg>/queries/],
    "docs":    [.md files in pkg],
    "config":  [paths under connections/ + environments/]
  },
  "features": [
    {
      "id":          "feature:<Pkg>:<slug>",
      "name":        "Human-readable feature name",
      "description": "What it does for the user",
      "members":     [ {kind, id, name, friendly_name, role, ...}, ... ],
      "ast_baseline": {
        "implements": [...repo-relative paths...],
        "tests":      [...repo-relative paths...]
      }
    }
  ]
}
```

## What to do

For each feature, decide which **files** in the package implement it and
which **test files** test it. Use the `members` list (it tells you the
RegisteredFunctions, Scripts, etc. that comprise the feature) plus the
file listings to find the real implementation. You may use `Read`, `Glob`,
and `Grep` tools to inspect any file under the package directory.

## Rules

- Only emit paths that appear in the `files` listing — exact match.
- Skip `package.g.ts` and any `*-api.ts` / `*.api.g.ts` (auto-generated).
- `package.ts` IS valid — it's the entry point that wires everything up.
- A feature typically has 2–10 implementation files and 0–3 test files.
- Be conservative: if a file is referenced but is just a one-liner type
  re-export, skip it. If you're not sure, skip.
- Tests must be tests of *this* feature (file mentions the feature's name,
  one of its functions, a related class, etc.) — not just any test in the
  package.
- For demo files (`src/demo*/`, `src/apps/`-like), if they specifically
  demo this feature, include them in `implements` with `role: "demo"`.

## Your output (STRICT JSON, no other text)

```json
{
  "package_id": "pkg:<Pkg>",
  "features": [
    {
      "id": "feature:<Pkg>:<slug>",
      "implements": [
        {"path": "public/packages/<Pkg>/src/widgets/foo.ts", "role": "core",  "confidence": 0.95},
        {"path": "public/packages/<Pkg>/src/utils/bar.ts",   "role": "support", "confidence": 0.85}
      ],
      "tests": [
        {"path": "public/packages/<Pkg>/src/tests/foo-tests.ts", "confidence": 0.9}
      ],
      "notes": "1-sentence justification, optional"
    }
  ]
}
```

`role` ∈ `{"core", "support", "demo", "config"}`. `confidence` ∈ [0, 1].

Return ONLY the JSON object. No markdown fences, no commentary.
