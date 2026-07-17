# Grokky changelog

## v.next

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