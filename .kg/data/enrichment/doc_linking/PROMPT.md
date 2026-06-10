# Doc → code linking for one help/ subtree

You are linking documentation pages in `public/help/` to the code entities
they describe. Your output drives `Documents` edges in the Datagrok
knowledge graph.

## Input

A JSON object with:
- `slice_id` — short id for this batch (e.g. `help-visualize`).
- `docs` — DocPage entries (`id`, `url`, `title`, `snippet` of body).
- `candidate_packages` — every Package in the corpus.
- `candidate_features` — every Feature (LLM-clustered groups of functions/scripts/docs that already hang together).
- `candidate_functions` — user-facing role-tagged Functions (role in {app, viewer, widget, panel, fileViewer, dashboard, semTypeDetector, transform, editor, …}).

Use IDs verbatim. Do NOT invent ids.

## Your output (STRICT JSON, no markdown fences, no commentary)

```json
{
  "slice_id": "help-visualize",
  "links": [
    {
      "doc_id": "doc:/help/visualize/viewers/scatter-plot",
      "target_id": "func:Charts:Sankey",         // or feature:..., or pkg:...
      "confidence": 0.92,
      "evidence": "title and content describe Sankey diagrams"
    }
  ]
}
```

## Rules

1. **Prefer Feature → DocPage links over Function → DocPage** when a feature
   exists. Features already group their member functions, so the indirect
   chain `Function ←PART_OF_FEATURE— Feature ←DOCUMENTS— DocPage` is more
   meaningful than spamming direct edges.

2. **Emit Function-level links only when the doc is unambiguously about a
   single function** (e.g. `/help/visualize/viewers/scatter-plot` → the
   scatter-plot viewer function specifically).

3. **Package-level links are useful for** package overviews, READMEs, or
   getting-started pages that cover many features at once
   (e.g. `/help/datagrok/solutions/domains/chem/overview` → `pkg:Chem`).

4. **A doc may link to multiple targets.** Common case: a "Cheminformatics
   solutions" doc links to `pkg:Chem` AND several `feature:Chem:*` features.
   Emit one row per target.

5. **Confidence**:
   - `0.95+` — title or URL slug exactly matches a Feature/Function name.
   - `0.80-0.94` — strong textual match (snippet describes the entity unambiguously).
   - `0.60-0.79` — domain match without specific naming alignment.
   - Below `0.60` — don't emit. Skip.

6. **Cross-domain links are fine.** A "Add new column" doc in `/help/transform/`
   may legitimately reference a Power Pack feature. Cross-package linking is
   the whole point of this enricher.

7. **One doc → many targets is normal.** Don't try to make it 1:1.

8. **Be exhaustive but not noisy.** Aim to link every doc to ≥ 1 target if
   any plausible match exists, but don't reach. Better to skip an unclear
   doc than to invent a wrong link.

9. **Reserve viewer-class slug matching** (e.g. `scatter-plot.md` → scatter-plot
   viewer). Datagrok core viewers (Bar Chart, Scatter Plot, Box Plot, …) are
   built into the platform itself — they may NOT have a corresponding
   RegisteredFunction in the input. In that case, link to the most relevant
   Charts/PowerGrid feature instead, with a slightly lower confidence (0.7-0.8).

Return ONLY the JSON. No prose around it.
