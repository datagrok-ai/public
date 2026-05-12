# Grokky changelog

## v.next

* AI: Fixed Enter key not submitting prompts in Chrome ≤ 50 / Dartium (no `KeyboardEvent.key`)
* AI: Prompts claimed by Datagrok's built-in handler now show a green "Handled natively" check instead of a Responses block
* AI: Ctrl+[ / Ctrl+] cycle through the current session's prompt history
* AI: Replaced the "Responses" accordion with a hover-only minimize icon that collapses a reply block to one line
* AI: Chat panel CSS now works in Chrome 50 — replaced flexbox `gap` with sibling margins, added `-webkit-` fallbacks for `filter` / `user-select` / `mask-image`
* AI: Unified the bottom controls row icons — same size, larger spacing, subtle hover background
* AI: Copy / thumbs-up / thumbs-down icons under a reply now appear only on hover of that message
* Added `src/polyfills.ts` (Chrome 50 / Dartium) — `crypto.randomUUID`, `Object.values`/`entries`/`fromEntries`, `String.prototype.trimStart`/`trimEnd`, `Array.prototype.flatMap`, `Element.prototype.append`/`prepend`/`replaceWith`; routed clipboard writes through a `copyToClipboard()` helper with an `execCommand` fallback

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