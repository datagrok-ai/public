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

1. **Author the handler as a decorated `static` method on `PackageFunctions` in `src/package.ts`** (`DG-FACT-179`; do not hand-edit `package.g.ts`). Signature must be `async (scriptCall: DG.FuncCall): Promise<void>` — use parameter name `scriptCall` (`DG-FACT-174`).

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
   Expected: `npx grok check` exits `0`; `package.g.ts` carries `//meta.role: scriptHandler` above the codegen-emitted function.

2. **Include all four mandatory annotations** — registration silently fails if any are missing (`DG-FACT-173`, `DG-FACT-175`): `meta.role: scriptHandler`, non-empty `meta.scriptHandler.language`, comma-separated `meta.scriptHandler.extensions` (no leading dot), and the single `funccall scriptCall` input.

3. **Add optional annotations as needed** — full list (`templateScript`, `codeEditorMode`, `commentStart`, `friendlyName`, `vectorizationFunction`, `parserFunction`, `meta.icon`) in `DG-FACT-176` and `DG-FACT-459`.

4. **If you set `vectorizationFunction`, register the target** with `//input: script script` + `//output: string result` — without these, lookup fails silently (`DG-FACT-177`).

   ```typescript
   // package.g.ts (emitted)
   //input: script script
   //output: string result
   export function makeVectorCode(script: any): string {
     return PackageFunctions.makeVectorCode(script);
   }
   ```

5. **Implement the handler body** — read `scriptCall.func.script`, execute, write outputs via `scriptCall.setParamValue` (or `scriptCall.outputs[...]`). The script's `#name/#input/#output` header is the source of truth for parameters (`DG-FACT-178`).

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

- **Function compiles but never registers** — a mandatory annotation is missing/malformed (case, leading dots, empty language). Startup skips silently (`DG-FACT-175`).
- **Hand-edited `package.g.ts` is overwritten** — edit the decorator in `package.ts` (`DG-FACT-179`).
- **Vectorization silently does nothing** — malformed nqName or missing `//input: script script` + `//output: string result` (`DG-FACT-177`).
- **Outputs blank** — handler returned a value instead of calling `scriptCall.setParamValue` (`DG-FACT-174`).
- **Wrong handler fires** — dispatch is by `#language:` header in the script body, not extension (`DG-FACT-178`).

## See also

- Source: `help/develop/how-to/scripts/custom-script-handlers.md`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-173` … `DG-FACT-179`, `DG-FACT-459` (annotations, dispatch,
  vectorization, decorator authoring); `DG-FACT-DRIFT-069`.
- Reference: `packages/Pyodide/src/package.ts:307-370` (decorator
  handler + vectorization function — only public-repo package using
  `scriptHandler`); `packages/Pyodide/src/package.g.ts:9-26` (codegen).
- Related skills: `python-functions`, `docker-containers`, `publish-packages`.
