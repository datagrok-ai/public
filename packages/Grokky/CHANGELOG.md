# Grokky changelog

## v.next

* AI: Fixed markdown tables rendering as glued pipe-text in chat — models emit `**Section**` directly followed by `| header |` rows, which the platform's spec-strict GFM parser refuses (a table cannot interrupt a paragraph); the renderer now inserts the missing blank line before table blocks (fence-aware), in both streaming and final render. Also added a paragraph break between pre-tool and post-tool text segments ("…check those.No, there…")
* AI: Gate blocks (grounding / verification) now run a background revision instead of double-answering or flashing drafts — the streamed answer stays visible, the revision streams hidden (runtime suppresses its chunks after the block), a "Revising…" status appears only if it takes a while, and the model decides the outcome: it writes a complete replacement (swapped in atomically on `final.revision: 'replaced'`) or replies `NO_REVISION` to keep the original untouched ('kept'). Block-reasons are reframed as internal feedback (never mention/quote it — fixes the leaked "the hook feedback doesn't apply…"), backed by an "Internal feedback" system-prompt rule
* AI: Help grounding is index-first — the runtime generates `workspace/help/INDEX.md` (path + title + description + keywords per published page, from the docs frontmatter; regenerated at startup and after each workspace sync), and the system prompt + GroundingGate now prescribe Read-the-index → Read-the-page instead of grepping 500+ files; WebFetch/WebSearch and `js-api/` reads now count as grounded, and the repo's `CLAUDE.md`/`.claude`/`.kg` are excluded from staged user workspaces (a measured help-turn timeout mechanism). Targets the measured 8–33 tool round-trips per help turn (see docs/LATENCY.md)
* AI: `datagrok_exec` self-verifies — a new optional `verify: {assertion, description}` parameter runs the verification assertion in the browser right after the action code (fresh scope, same round-trip) and returns it as `verified`; a passing one satisfies the Verifier, dropping the separate `datagrok_verify` model round-trip that made 2 round-trips the floor for every action turn
* AI: Enriched all ~440 published help pages with frontmatter `keywords` (searchable synonyms and task phrases) and one-line `description`s; ~110 legacy pages got proper frontmatter with their body H1 promoted to `title`
* AI: Added a latency/accuracy benchmark harness (`docs/BENCHMARK.md`) to put numbers on every LATENCY.md lever — `Grokky:runBenchmark(label)` drives a golden prompt suite (`files/benchmark/suite.yaml`, ~20 prompts across help/visualization/analysis/codegen/multitool/query) through the real runtime pipeline, timing TTFT/total/tool-round-trips per turn, reading the SDK token/cost usage now forwarded on the `final` event, scoring each prompt (deterministic assert + Haiku judge), and downloading a JSON+Markdown report; `Grokky:compareBenchmarks(a, b)` diffs two runs into a delta table. Latency/context metrics work immediately; token/cost columns need a claude-runtime image rebuild
* AI: Prompt suggestions (empty-state cards + wand icon) now appear on any view, not only table views — the always-applicable "Anywhere" and "Code generation" blocks show from the home screen; view/column blocks still self-add when a table is open
* AI: Expanded and diversified the table-view prompt suggestions (empty state + wand menu) — added semantic blocks for sequences, 3D structures, potency (IC50/EC50/concentration), geospatial, money and free text, plus themed whole-table blocks (Visualize, Analyze & model, Time trends, Clean & transform); blocks and individual suggestions self-hide when the data lacks the required column, and a new `isDateTime` slot filter enables time-series suggestions
* AI: Ambiguous suggestions no longer make the assistant guess — parameter-dependent ones (substructure/motif search, similarity reference, calculated-column formula, custom viewer, choropleth metric, data search) now post their clarifying question straight into the chat as the assistant's reply via a new `immediateResponse` suggestion field (no AI round-trip); the user's answer then runs the task with the intent carried forward as context. "Open a demo dataset" uses a new client-side `action` field to render an inline choice block in the panel (Demographics / SPGI / curves / beer / random walk) that reads like an assistant reply; picking one removes the block, opens the table, and posts an "Opened …" note, and the block is dismissed if the user moves on to a new prompt
* AI: Views brief the assistant — new `View.aiDescription` (js-api `ViewBase`/`View`, Dart `View` with JsViewHost forwarding) carries a short what-this-view-is + which-functions-to-call note, prepended to every prompt as "About this view"; query and script editors ship one
* AI: `list_view_functions` matching is now OR-ranked (functions matching more query words rank higher) — AND-matching returned nothing for multi-word queries; a zero-match response now says how many functions the view actually has
* AI: View functions replace view AI tools — the assistant now works with the platform's `getFunctions()` (registered `DG.Func`s applicable to the view: `meta.viewType` functions, Dart view overrides, JS views via `view.jsView`) through three static meta-tools (`list_view_functions` search with ≤10 results, `get_view_function_result`, `call_view_function`), so views with hundreds of functions no longer bloat the context; `getAITools`/`AIViewTool`/`viewAIToolsProvider` are removed
* AI: Restored empty-state prompt suggestions on table views — the singleton panel now resolves suggestion cards against the live current view (was pinned to its creation view) and refreshes them on view switches
* AI: Query editor — removed the dedicated DB panel and catalog selector; the assistant now works on the query view through tools: `get_query_info` / `set_query_and_run` (Dart-native) plus SQL schema exploration and test-execution tools (`list_db_*`, `get_db_table_details`, `get_sql_test_result`)
* AI: Script editor — removed the dedicated scripting panel and language selector; the assistant reads/writes the open editor via `get_script_code` / `set_script_code` and infers the language from the script header

* AI: The AI panel is now a single persistent assistant — switching views no longer swaps panels or resets the conversation; workspace context (current view, all open views, all tables) is rebuilt fresh on every prompt
* AI: `datagrok-exec` / `datagrok_verify` blocks run against the live current view, so code targets views Claude just opened (e.g. a joined table) instead of the view the prompt started from
* AI: Query-editor AI assistant no longer force-opens the AI panel when a query view opens — it registers and shows only on its toggle icon; it's disposed with its view
* AI: Script-generation panel became a singleton rebound to the script view whose AI icon invoked it — one scripting conversation across script views

* GROK-18695: Dependency security updates — refreshed lockfile so dev-only puppeteer no longer pins the runtime ws at a vulnerable 8.17.0 (DoS/memory-disclosure advisories)
* GROK-18695: Docker: raised security floors (VEX) — claude-runtime: deb12 upgrade, npm@latest refresh, pip/setuptools bump, hono/@hono/node-server/ws floors, dropped unused workspace JDBC jars; mcp-server: express-rate-limit/fast-uri/path-to-regexp/qs/ip-address overrides
* GROK-20054: Report: Error: Claude runtime container is not running
* AI: Fixed query/script view going blank (tab unselected) on first click when switching back to it with the AI panel open — defer the panel dock/undock until after the view switch settles
* AI: Replaced the 15-min skill/agent sync poll with an on-demand, TTL-gated refresh of packages + shared connections — idle sessions no longer sync
* AI: Fixed shared-connection files that fail to download (e.g. invalid JSON) re-downloading on every sync
* AI: Voice input — while a prompt is being processed the loader shows a `Say "cancel" to stop` hint, and saying "stop"/"cancel" aborts the run instead of being sent as a new prompt
* AI: Fixed Enter key not submitting prompts in Chrome ≤ 50 / Dartium (no `KeyboardEvent.key`)
* AI: Prompts claimed by Datagrok's built-in handler now show a green "Handled natively" check instead of a Responses block
* AI: Ctrl+[ / Ctrl+] cycle through the current session's prompt history
* AI: Replaced the "Responses" accordion with a hover-only minimize icon that collapses a reply block to one line
* AI: Chat panel CSS now works in Chrome 50 — replaced flexbox `gap` with sibling margins, added `-webkit-` fallbacks for `filter` / `user-select` / `mask-image`
* AI: Unified the bottom controls row icons — same size, larger spacing, subtle hover background
* AI: Copy / thumbs-up / thumbs-down icons under a reply now appear only on hover of that message
* Added `src/polyfills.ts` (Chrome 50 / Dartium) — `crypto.randomUUID`, `Object.values`/`entries`/`fromEntries`, `String.prototype.trimStart`/`trimEnd`, `Array.prototype.flatMap`, `Element.prototype.append`/`prepend`/`replaceWith`; routed clipboard writes through a `copyToClipboard()` helper with an `execCommand` fallback
* AI: Added a Run button to the ribbon of file views opened from `MyFiles/agents/scripts/`
* AI: Added wand icon to open curated prompt suggestions menu (loaded from suggestions.yaml)

## 1.0.4 (2026-02-01)

* Grokky: Azure compatible client
* AI: Made a namespace, made a nested config, improved namings
* Got rid of "latest" versions
* Grokky: Refactor to Vercel AI

## 1.0.1 (2025-12-18)

* Grokky: Move stuff to proxy

## 1.0.0 (2025-12-17)

### Features

* Grokky: Invoking function and chains of functions
* Grokky: Fuzzy matching
* AI: Powersearch Search AI features
* Grokky: Add structured output and reuse Gemini session to avoid recreating on each call
* Grokky: AI sql

### Bug fixes

* GROK-19198: Duplicating AI icons issue