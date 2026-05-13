---
name: custom-script-handlers
version: 0.1.0
description: |
  Wire a Datagrok package function as the executor for a scripting
  language the platform doesn't ship natively (Clojure, an internal
  DSL, a WASM runtime) so the language appears in the "Create a
  script" picker and every script whose `#language:` header matches
  your tag routes to your code. For plugin authors plugging a new
  language or sandboxed runtime into the platform without forking it.
  Produces a decorated `async (scriptCall: DG.FuncCall): Promise<void>`
  method on `PackageFunctions` plus the matching `package.g.ts` codegen.
  Use when asked to "plug a new language into the platform", "run a
  sandboxed interpreter for scripts users author", or "make
  `#language: <ours>` scripts dispatch to my package".
triggers:
  - plug a new language into datagrok
  - register an interpreter for user scripts
  - route a custom dsl to my package
  - add a runtime for an internal dsl
  - run wasm code from a script header
  - dispatch scripts by language tag
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# custom-script-handlers

## When to use

Your package needs Datagrok to treat a not-yet-supported scripting
language as first-class — show it in "Create a script", route every
script whose `#language:` header matches your tag to your code, and
optionally auto-vectorize. For Python / R / JS / Julia, use the
built-in runtimes — no handler needed.

## Prerequisites

- Decorator-based package (`grok create <Name>`): `PackageFunctions`
  class in `src/package.ts`, codegen emits `src/package.g.ts`.
- An execution mechanism reachable from JS — in-process interpreter,
  WebAssembly module, or Docker container (`docker-containers` skill).
- Familiarity with `DG.FuncCall` (`scriptCall.inputs`, `.outputs`,
  `.setParamValue`, `.func.script`).

## Steps

1. **Author the handler as a decorated `static` method on `PackageFunctions` in `src/package.ts`.**
   Do NOT hand-edit `package.g.ts` — codegen overwrites it every
   build (`DG-FACT-179`). Signature MUST be
   `async (scriptCall: DG.FuncCall): Promise<void>`; name the
   parameter `scriptCall` (not `call`) — every shipping handler and
   the codegen emit `scriptCall` (`DG-FACT-174`).

   ```typescript
   // src/package.ts
   export class PackageFunctions {
     @grok.decorators.func({
       meta: {
         'role': 'scriptHandler',
         'scriptHandler.language': 'clojure',
         'scriptHandler.extensions': 'clj,cljs,cljr',
         'scriptHandler.commentStart': ';',
         'scriptHandler.codeEditorMode': 'clojure',
         'icon': 'files/clojure.png',
       },
     })
     static async clojureScriptHandler(scriptCall: DG.FuncCall): Promise<void> {
       // 1. read body:     scriptCall.func.script
       // 2. execute (in-process / WASM / Docker)
       // 3. write outputs: scriptCall.setParamValue('<name>', value);
     }
   }
   ```
   Expected: `npx grok check` exits `0`; `src/package.g.ts` contains a
   matching `//meta.role: scriptHandler` block above an
   `async function clojureScriptHandler(scriptCall: DG.FuncCall): Promise<void>`
   (shape per `packages/Pyodide/src/package.g.ts:15-26`).

2. **Include all four mandatory annotations — registration silently fails if any are missing.**
   Startup only registers a function when `meta.role` equals literal
   `scriptHandler` AND `meta.scriptHandler.language` is non-empty
   (`DG-FACT-173`, `DG-FACT-175`):
   - `meta.role: scriptHandler` (case-sensitive;
     `DG.FUNC_TYPES.SCRIPT_HANDLER`).
   - `meta.scriptHandler.language: <name>` — what scripts put in their
     `#language:` header to dispatch here (`DG-FACT-178`).
   - `meta.scriptHandler.extensions: <ext1,ext2>` — comma-separated,
     **no leading dot** (`py`, not `.py`); split on `,` at startup.
   - One input `funccall scriptCall` (the only parameter).

3. **Add optional annotations as needed — seven exist** (`DG-FACT-176`).
   - `scriptHandler.templateScript` — boilerplate for new scripts;
     embed newlines as literal `\n`.
   - `scriptHandler.codeEditorMode` — CodeMirror-5 mode (`python`,
     `clojure`); defaults to `language`.
   - `scriptHandler.commentStart` — line-comment chars; default `#`.
   - `scriptHandler.friendlyName` — UI label; defaults to `language`.
   - `scriptHandler.vectorizationFunction` — nqName `<Pkg>:<fn>` of a
     `DG.Script → string` function (step 4).
   - `scriptHandler.parserFunction` — nqName of a text→`DG.Script`
     parser (`DG-FACT-459`); rarely needed (built-in `#name:` parser
     covers Python / Clojure / JS-style headers).
   - `meta.icon` — package-relative path (e.g. `files/clojure.png`).

4. **If you set `vectorizationFunction`, register the target with the right input/output types.**
   Target MUST be a package function with signature
   `(script: DG.Script) => string`. Codegen needs `//input: script script`
   + `//output: string result`; without these the lookup fails silently
   at handler-construction (`DG-FACT-177`). Reference
   (`packages/Pyodide/src/package.g.ts:9-13`):

   ```typescript
   // package.g.ts (emitted)
   //input: script script
   //output: string result
   export function makeVectorCode(script: any): string {
     return PackageFunctions.makeVectorCode(script);
   }
   ```

5. **Implement the handler body — read script, execute, write outputs via `scriptCall`.**
   The script body is on `scriptCall.func.script` (a `DG.Script`); its
   declared inputs/outputs come from the script's own
   `#name/#input/#output` header, NOT the handler's annotations
   (`DG-FACT-178`). Read inputs as `scriptCall.inputs[<name>]`; write
   outputs via `scriptCall.setParamValue('<name>', value)` (or
   `scriptCall.outputs[<name>] = value`). Reference
   (`packages/Pyodide/src/package.ts:361-369`):

   ```typescript
   static async pyodideLanguageHandler(scriptCall: DG.FuncCall): Promise<void> {
     const req = await prepareRequest(scriptCall);   // reads .func.script + .inputs
     const response = await sendRequest(req);
     await setOutputs(scriptCall, response);         // .setParamValue(...)
   }
   ```

6. **Build, regenerate codegen, publish.**
   ```bash
   npm install
   npx grok check                  # exits 0; package.g.ts regenerated
   npx grok publish <host>         # add --release once stable
   ```
   Expected: **Tools → Scripting → New script** shows your
   `friendlyName` (or `language`); scripts whose header is
   `#language: <yourLanguage>` dispatch to your handler.

## Common failure modes

- **Function compiles but never registers.** Mandatory annotation
  missing or malformed — `meta.role` typo (case-sensitive
  `scriptHandler`), extensions with leading dots (`.py` vs `py`), or
  empty `scriptHandler.language`. Startup silently skips it
  (`DG-FACT-175`); `npx grok check` still passes. Inspect
  `package.g.ts` — all four mandatory lines must be present.
- **Hand-edited `package.g.ts` is overwritten on next build.** Edit
  the decorator in `package.ts`; codegen regenerates `package.g.ts`
  from it every build (`DG-FACT-179`).
- **Vectorization silently does nothing.** Either the
  `vectorizationFunction` nqName is malformed (`Package:funcName`,
  case-sensitive) or the target lacks `//input: script script` /
  `//output: string result` and codegen never registered it
  (`DG-FACT-177`). No build error.
- **Script runs but outputs are blank.** Handler didn't write back.
  Outputs go through `scriptCall.setParamValue('<name>', value)` (or
  `scriptCall.outputs[<name>] = value`); returning a value does
  nothing — return type is `Promise<void>` (`DG-FACT-174`).
- **Script with the right extension still falls back to a built-in handler.**
  Dispatch is by `#language:` header in the script body, NOT by file
  extension (`DG-FACT-178`). Confirm the script header.

## See also

- Source: `help/develop/how-to/scripts/custom-script-handlers.md`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-173` … `DG-FACT-179`, `DG-FACT-459` (annotations, dispatch,
  vectorization, decorator authoring); `DG-FACT-DRIFT-069`.
- Reference: `packages/Pyodide/src/package.ts:307-370` (decorator
  handler + vectorization function — only public-repo package using
  `scriptHandler`); `packages/Pyodide/src/package.g.ts:9-26` (codegen).
- Related skills: `python-functions`, `docker-containers`, `publish-packages`.
