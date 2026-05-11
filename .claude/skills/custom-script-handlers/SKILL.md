---
name: custom-script-handlers
description: Register a Datagrok package function as a custom script-language handler so the platform recognises and executes scripts written in that language
---

# custom-script-handlers

## When to use

Your package needs Datagrok to treat a not-yet-supported scripting
language as first-class — show it in the "Create a script" dropdown,
route any script whose `#language:` header matches your tag to your
code, and (optionally) auto-vectorize it. Triggers: "add Clojure
support", "run Pyodide-style Python in the browser", "register a
custom script runtime". For Python/R/JS scripts on the built-in
runtimes, no handler is needed — author the script directly.

## Prerequisites

- A package scaffold (`grok create <Name>`) using decorator-based
  authoring (`PackageFunctions` static class in `src/package.ts`,
  codegen emits `src/package.g.ts`).
- An execution mechanism reachable from JS — in-process interpreter,
  WebAssembly module, or Docker container (`docker-containers`).
- Familiarity with `DG.FuncCall` (`scriptCall.inputs`,
  `scriptCall.outputs`, `scriptCall.setParamValue`,
  `scriptCall.func.script`).

## Steps

1. **Author the handler as a decorated `static` method on `PackageFunctions` in `src/package.ts`.**
   Do NOT hand-write the comment-header form in `package.g.ts` —
   codegen overwrites it (`DG-FACT-179`). The signature MUST be
   `async (scriptCall: DG.FuncCall): Promise<void>` (`DG-FACT-174`);
   use the name `scriptCall`, not `call`, to match codegen and every
   shipping handler (`DG-FACT-DRIFT-070`).

   ```typescript
   // src/package.ts
   import * as DG from 'datagrok-api/dg';
   import * as grok from 'datagrok-api/grok';
   export * from './package.g';
   export const _package = new DG.Package();

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
       // 1. read body: scriptCall.func.script
       // 2. execute (in-process / WASM / Docker)
       // 3. write each output: scriptCall.setParamValue('<name>', value);
     }
   }
   ```
   Expected: `npx grok check` exits `0`; `src/package.g.ts` contains a
   matching `//meta.role: scriptHandler` header above an
   `async function clojureScriptHandler(scriptCall: DG.FuncCall): Promise<void>`.

2. **Include all four mandatory annotations — registration silently fails if any are missing.**
   Startup only registers a function when `meta.role` equals the
   literal `scriptHandler` AND `meta.scriptHandler.language` is
   non-empty (`DG-FACT-173`, `DG-FACT-175`):
   - `meta.role: scriptHandler` (literal token).
   - `meta.scriptHandler.language: <name>` — what scripts put in their
     `#language:` header to dispatch here (`DG-FACT-178`).
   - `meta.scriptHandler.extensions: <ext1,ext2,...>` — comma-separated,
     **no leading dot** (`py`, not `.py`); split on `,` at startup.
   - One input `funccall scriptCall` (the only parameter).

3. **Add optional annotations as needed — seven exist; the article documents only four** (`DG-FACT-176`, `DG-FACT-DRIFT-069`).
   - `scriptHandler.templateScript` — boilerplate inserted into newly
     created scripts. Embed newlines as literal `\n`.
   - `scriptHandler.codeEditorMode` — CodeMirror-5 mode name
     (`python`, `clojure`); defaults to `language`.
   - `scriptHandler.commentStart` — line-comment character; default
     `#`. **Article-undocumented; used by Pyodide.**
   - `scriptHandler.vectorizationFunction` — nqName `<Pkg>:<fn>` of a
     `DG.Script → string` function (see step 4).
   - `scriptHandler.friendlyName` — UI label override; defaults to
     `language`. **Article-undocumented.**
   - `scriptHandler.parserFunction` — nqName of a text→`DG.Script`
     parser. **Article-undocumented.**
   - `meta.icon` — package-relative path (e.g. `files/clojure.png`).

4. **If you set `vectorizationFunction`, register the target function with the correct input/output types.**
   The targeted function MUST be a package function in the same
   package, signature `(script: DG.Script) => string`. Codegen needs
   `//input: script script` + `//output: string result` (`DG-FACT-177`,
   `DG-FACT-DRIFT-071`); without these the lookup fails silently at
   handler-construction. Reference (`packages/Pyodide/src/package.g.ts:9-13`):

   ```typescript
   // package.g.ts (emitted)
   //input: script script
   //output: string result
   export function makeVectorCode(script: any): string {
     return PackageFunctions.makeVectorCode(script);
   }
   ```

5. **Implement the handler body — read the script, execute it, write outputs back via `scriptCall`.**
   The script body is on `scriptCall.func.script` (a `DG.Script`) and
   its declared inputs/outputs come from the script's own
   `#name/#input/#output` header, NOT from the handler's annotations
   (`DG-FACT-178`). Pass each input via `scriptCall.inputs[<name>]`;
   write each output via `scriptCall.setParamValue('<name>', value)`
   (or assign to `scriptCall.outputs[<name>]`). Reference
   (`packages/Pyodide/src/package.ts:361-369`):

   ```typescript
   static async pyodideLanguageHandler(scriptCall: DG.FuncCall): Promise<void> {
     const req = await prepareRequest(scriptCall);   // reads scriptCall.func.script + .inputs
     const response = await sendRequest(req);        // executes
     await setOutputs(scriptCall, response);         // scriptCall.setParamValue(...)
   }
   ```

6. **Build, regenerate codegen, publish.**
   ```bash
   npm install
   npx grok check                  # exits 0; package.g.ts regenerated
   npx grok publish <host>         # add --release once stable
   ```
   Expected: after publish, "Create a script" → language dropdown
   shows your `friendlyName` (or `language`); any script whose header
   contains `#language: <yourLanguage>` dispatches to your handler.

## Common failure modes

- **Function compiles but never registers.** A mandatory annotation
  is missing or malformed — `meta.role` typo (case-sensitive
  `scriptHandler`), extensions with leading dots (`.py` instead of
  `py`), or `scriptHandler.language` empty. Startup silently skips
  it (`DG-FACT-175`). `npx grok check` still passes, but the
  language is absent from the dropdown. Inspect the generated
  `package.g.ts` block — all four lines must be present.
- **Hand-edited `package.g.ts` is overwritten on next build.** Edits
  belong in `package.ts` decorators; codegen regenerates
  `package.g.ts` from them every build (`DG-FACT-179`).
- **Vectorization silently does nothing.** Either `vectorizationFunction`
  nqName is malformed (`Package:funcName`, case-sensitive) or the
  target function is missing `//input: script script` / `//output:
  string result` and codegen never registered it (`DG-FACT-177`,
  `DG-FACT-DRIFT-071`). No build error; lookup fails at
  handler-construction.
- **Script runs but outputs are blank.** The handler didn't write
  back. Outputs go through `scriptCall.setParamValue('<name>', value)`
  (or `scriptCall.outputs[<name>] = value`); returning a value does
  nothing — return type is `Promise<void>` (`DG-FACT-174`).
- **Script with the right extension still falls back to the built-in handler.**
  Dispatch is by `#language:` header inside the script body, NOT by
  file extension (`DG-FACT-178`). Check the script's first lines
  include `#language: <yourLanguage>`.

## Verification

- `npx grok check` exits `0`; `src/package.g.ts` contains all four
  mandatory annotations above your handler function.
- After `grok publish`, **Tools → Scripting → New script** dropdown
  lists your language; selecting it loads the `templateScript` body
  (if set) and switches the editor to `codeEditorMode`.
- Run a minimal script (`#language: <yours>` + one `#input` + one
  `#output`) from the editor; the output panel shows the value set
  via `scriptCall.setParamValue`.
- If `vectorizationFunction` is set: tabular invocation (apply to a
  column) does NOT fall back to per-row scalar invocation.
- Browser console — no `scriptHandler` warnings or class-load errors.

## See also

- Source articles:
  - `help/develop/how-to/scripts/custom-script-handlers.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-173` … `DG-FACT-179` (role token, mandatory + optional
  annotations, vectorization signature, dispatch rule, decorator
  authoring) and drifts `DG-FACT-DRIFT-069` … `DG-FACT-DRIFT-071`
  (article-undocumented annotations; canonical `scriptCall` param
  name; vectorization target's required annotations).
- Reference packages:
  - `packages/Pyodide/src/package.ts:307-370` — decorator-form
    handler + vectorization function; the only public-repo package
    implementing `scriptHandler`.
  - `packages/Pyodide/src/package.g.ts:9-26` — codegen output to
    expect after `grok check`.
- Related skills:
  - `python-functions` — built-in runtime; no handler needed.
  - `docker-containers` — execution mechanism for handlers that
    shell out.
  - `publish-packages` — final `grok publish` step.
