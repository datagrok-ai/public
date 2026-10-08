# CLAUDE.md - Reactive Tree Driver

**Important**: This library is the core engine for `Compute2`. For full development context, build instructions, testing, and publishing guidelines, see [`packages/Compute2/CLAUDE.md`](../../../packages/Compute2/CLAUDE.md).

## Quick Reference

- **RTD source**: `libraries/compute-utils/reactive-tree-driver/`
- **UI consumer**: `packages/Compute2/` (Vue 3 components)
- **Tests**: `packages/LibTests/src/tests/compute-utils/reactive-tree-driver/`
- **Build**: `cd libraries/compute-utils && npm run build` (runs `tsc`, outputs `.js` + `.d.ts` alongside `.ts` files)
- **Build everything**: `cd packages/Compute2 && npm run build-all`
- **Publish**: Always use `grok publish --release` for compute packages
- **User docs**: `help/compute/workflows/` — configuration, link types, link spec language, code usage, examples

## What This Library Does

Reactive tree driver propagates data through dynamically created and mutated function call trees. It manages:
- Tree state (static/dynamic pipeline nodes, action steps, FuncCall leaves)
- Data/validator/meta links between steps
- Consistency tracking and validation
- Serialization/deserialization of pipeline state

### Action Steps

`type: 'action'` is a config-level step type that gets converted to a `PipelineConfigurationStaticProcessed`
with `isActionStep: true` during config processing. It has no children, no links, no history. Designed as
a `visibleOn` target for actions from outer pipelines. It may carry an optional `nqName` pointing at a
function whose `help`/`readme` option supplies its help panel content (it has no FuncCall of its own).
The type guard `isPipelineActionConfig()` identifies action configs before processing; after processing
they are regular static pipeline configs.

## Rules, Checks and Annotations Internals

User docs (`help/compute/workflows/rules-and-checks.mdx`) describe behaviour only; the mechanics live here.

- **Formulas** (`src/config/rule-formula.ts`): `expandLinks` compiles formula strings into the existing JSON Logic, effect and source objects before `expandRule`/`expandCheck`, so expansion and the runtime only see objects. A string compiles in `when` (rule, effect, check), each `effects` element and each `sources` entry; in value fields (`items`, `message`, `value`, `values`, `meta.*`, object-source `args`, check `message`) only when it starts with `=`, else it stays literal; strings nested in objects never compile. Translation beyond a 1:1 rewrite: symbolic ops are renamed (`gt` → `>`), the alias arguments of `missing`/`missing_some`/`var` become path strings, `true`/`false`/`null` are keywords (an alias with such a name is read as `var("null")`), the element of `map`/`filter`/`all`/`some`/`none` is the driver name `$it` (`{var: ''}`, an error outside the element argument), `obj()` compiles to `{literal}` with constant values only (`evaluate` swaps each `{literal}` for a `$literal` op reading the active literal list, so it also works in element arguments, where JSON Logic replaces the data with the element), calls always emit argument arrays, expression ops take positional arguments only, and a double-quoted string accepts only the `\"`, `\\`, `\n`, `\t` escapes so a regex pattern is not silently stripped. `formulaOps` takes the driver's op names from `driverOpNames` (rule-expressions.ts), so a new driver op is callable from formulas at once. The effect and source call tables are typed over `RuleEffect['effect']` and the `RuleSource` kinds, so a new kind fails the build until it gets a call entry, and every field of an effect type must be listed in its call (`_EveryEffectFieldHasACall`); `js` and dataframe `table` sources stay objects. Malformed object shapes pass through untouched for `expandRule` to report. `expandRule` rejects a rule alias read as `{var: 'name'}` inside an element argument (`usedAliases` collects those into `fields`), where JSON Logic reads the element's field; `var("name")` (`{var: ['name']}`) is the explicit field read.
- **Rule expansion** (`src/config/rule-expansion.ts`): rule `r` becomes up to three links, `r::meta`, `r::validator`, `r::data`, one per effect family. Each copies `from`, `not`, `base`, `nodePriority`, `dataFrameMutations`, keeps only the `to` entries its effects target, and carries `when`/`effects`/`sources` in `params`. `debounce` goes to the validator link, `runOnInit` to the data link. `(call)` is not allowed in rule queries; the expansion adds the `$call` input that `validators`/`choices` sources without `names` need.
- **Sources** (`src/runtime/rule-sources.ts`): every expanded link resolves its sources on each run before `when`, so one rule may resolve a source up to three times. `file`, `table`, and `func`/`query` whose `args` read no alias (or that have none) are loaded once per link (`isConstant`). Synchronous results (no call made, sync `js`) keep the link batchable. `choices` returns `{items, values, inList, row, rowErrors}` via `FuncCall.evalParamChoices`, cached per call and io until a `dependsOn` input changes, a run that finds an evaluation in flight waits for it and checks its `dependsOn` then; `row` cells are converted to the input types like the form's editors parse text; datetime cells arrive as platform objects (js-api `toJs` does not convert nested lookup values) and go through `DG.toJs`, text and timestamps through `dayjs`. `$nonscalar` keeps datetime ios, unlike `DG.TYPES_SCALAR`.
- **Checks** (`src/config/checks.ts`): annotation options and `check` links share `expandChecks`; a `check` `c` becomes `c::min`, `c::table`, ... with `value:<io>` (and `$table:<query>`) inputs, plus the `vars` aliases (`<name>:<query>`) for expression options, and `$target:<io>` output; `visible` becomes a meta link; an `enabled:` annotation is folded into it (`(visible) && (enabled)`), so workflows hide what the form disables. Driver-added names carry a `$` prefix (`$call`, `$table`, `$target`, `$verdicts`, the `$all` expression key, the `$literal` op); the alias grammar accepts a leading `$`, and rules and checks reject it in user aliases, `vars` and source names. Only `value` stays shared with step inputs, matching the platform's validator expressions.
- **Annotation rules** (`annotationRules` in `rule-expansion.ts`): evaluated choices become `::<io>:choices` (items, and a `rowErrors` warning on the lookup key; no check of the value, as a `choices` annotation of any kind only fills the dropdown and a workflow may replace its items, so validation is an explicit `check` link or rule; no `emptyChoice` — the empty option follows the `nullable`/`optional` annotation in `InputForm`) and, for the first `propagateChoice: all` input only, `::<io>:lookup` (`runOnInit`, `assign` with `ignoreCase` and `restricted` into `inputs(nq, key|$nonscalar|$linked)`). They must stay expressible as hand-written rules.

## Testing Rules

Tests live in `packages/LibTests/src/tests/compute-utils/reactive-tree-driver/`, NOT alongside the RTD source.

**Mandatory conventions for all new RTD tests:**

- **RxJS virtual time only.** Use `TestScheduler` from `rxjs/testing` via `createTestScheduler()` (from `test-utils.ts`). Never use `async/await`, `setTimeout`, `toPromise()`, or real async operations.
- **Mock state only.** Use `mockMode: true` when creating `StateTree`. Use `FuncCallMockAdapter` / `MemoryStore` — never real `DG.FuncCall` instances or server calls.
- **Inline configs.** Define `PipelineConfiguration` objects directly in the test file. Do not fetch configs from the server via `callHandler()` or `grok.functions.call()`.
- **Marble diagrams for reactivity.** Test observable sequences with `expectObservable()`, `cold()`, `hot()` inside `testScheduler.run()`. Use standard marble notation: `'-a'`, `'a b'`, `'^ 1000ms !'`, etc.
- **Snapshot compare for tree structure.** Use `snapshotCompare(tree.toSerializedState({disableNodesUUID: true}), 'TestName')` for structural assertions.
- **`before()` runs once per category**, not per test. Reset any shared state inside each test or use `createTestScheduler()` which auto-resets frame/index.

See `packages/Compute2/CLAUDE.md` for full test categories and build/run instructions.

## After Making Changes

1. **Write tests for your changes.** Every new feature, bug fix, or behavior change must have corresponding tests. Add them to the appropriate file in `packages/LibTests/src/tests/compute-utils/reactive-tree-driver/`, or create a new file and register it in `packages/LibTests/src/package-test.ts`.
2. Rebuild compute-utils: `cd libraries/compute-utils && npm run build`
3. Rebuild LibTests: `cd packages/LibTests && npm run build-all`
4. Publish LibTests: `cd packages/LibTests && grok publish --release`
5. **Run ALL RTD tests** — not just the ones you added. Changes to the driver can break unrelated tests due to shared reactive state and link recalculation. Verify the full suite passes:
   ```bash
   grok test --skip-build --skip-publish --category "ComputeUtils: Driver"
   ```
   Environment-specific settings (Puppeteer path, `--host`) vary by machine — check Claude Code memory files for your local setup.
