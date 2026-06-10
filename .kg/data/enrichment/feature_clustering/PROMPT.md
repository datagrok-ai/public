# Feature clustering for one Datagrok plugin

You are analyzing a single Datagrok plugin (package). Your job is to identify
**Features** — coherent, user-facing capabilities that this plugin provides.

A Feature is NOT one function. It is a *theme* that may be implemented by:
- one or more `RegisteredFunction`s (the user-callable functions)
- supporting `Script`s (server-side compute)
- backing `DataQuery`s (SQL / data sources)
- a `Tutorial` that walks a user through it
- one or more `DocPage`s that explain it
- changelog entries describing its evolution

## Input

You will receive a JSON object describing the package: its metadata, its
functions (with role/tags/inputs/outputs/description/help_url/meta), scripts,
queries, tutorials, doc pages, properties, and a recent changelog. Every entity
has an `id`. Use those IDs verbatim in your output.

## Your output (STRICT JSON, no other text)

```json
{
  "package_id": "pkg:<Folder>",
  "features": [
    {
      "id": "feature:<Folder>:<slug>",
      "name": "Human-readable feature name",
      "description": "1-2 sentences describing what this feature gives the user.",
      "category": "<short tag, e.g. chem-search | chem-rendering | chem-mpo | viz | ml | data-import | infra>",
      "confidence": 0.0,
      "members": [
        {"id": "func:Pkg:foo",       "weight": 1.0,  "role": "core"},
        {"id": "script:Pkg:bar",     "weight": 0.7,  "role": "core"},
        {"id": "doc:/help/.../X",    "weight": 0.9,  "role": "doc"},
        {"id": "tutorial:Pkg:Foo",   "weight": 0.95, "role": "doc"}
      ],
      "evidence_summary": "1 sentence on why these belong together (shared name root, shared menu path, shared tutorial, etc.)"
    }
  ],
  "feature_relations": [
    { "from": "feature:Pkg:foo", "to": "feature:Pkg:bar", "kind": "subfeature|sibling|related" }
  ],
  "unassigned_member_ids": ["..."]
}
```

## Rules

- `id` slugs are **kebab-case**, derived from the feature name. Example: `"Activity Cliffs"` → `feature:Chem:activity-cliffs`.
- A feature may contain entities of any kind in `members`. Use only IDs from the input.
- `weight` ∈ [0, 1]. Use 1.0 for clearly core members, 0.6–0.9 for supporting members, lower for loose associations.
- `role` ∈ `{"core", "doc", "test", "sample", "related"}`.
- `confidence` ∈ [0, 1]. Use higher when multiple naming/tag/menu signals agree, lower when you're inferring.
- Be **conservative** about cross-feature relations. Emit `feature_relations` only when it's clearly subfeature/sibling/related.
- It is OK to have many features — typical plugin has 5–25 features. Don't lump unrelated functions together.
- Anything genuinely uncategorized goes in `unassigned_member_ids`. Be honest about uncertainty.
- Aim for features that a developer (or a user) would actually search for: "Substructure search", "MPO profiles", "Reaction enumeration", "Activity cliffs", "Scaffold tree", not "Chem widget #3".

Return ONLY the JSON object. No markdown fences, no commentary.
