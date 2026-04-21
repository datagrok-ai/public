# CLAUDE.md - Reactive Tree Driver

**Important**: This library is the core engine for `Compute2`. For full development context, build instructions, testing, and publishing guidelines, see [`packages/Compute2/CLAUDE.md`](../../../packages/Compute2/CLAUDE.md).

## Quick Reference

- **RTD source**: `libraries/compute-utils/reactive-tree-driver/`
- **UI consumer**: `packages/Compute2/` (Vue 3 components)
- **Tests**: `packages/LibTests/src/tests/compute-utils/reactive-tree-driver/`
- **Build**: `cd libraries/compute-utils && npm run build` (runs `tsc`, outputs `.js` + `.d.ts` alongside `.ts` files)
- **Build everything**: `cd packages/Compute2 && npm run build-all`
- **Publish**: Always use `grok publish --release` for compute packages

## What This Library Does

Reactive tree driver propagates data through dynamically created and mutated function call trees. It manages:
- Tree state (static/parallel/sequential pipeline nodes, action steps, FuncCall leaves)
- Data/validator/meta links between steps
- Consistency tracking and validation
- Serialization/deserialization of pipeline state

### Action Steps

`type: 'action'` is a config-level step type that gets converted to a `PipelineConfigurationStaticProcessed`
with `isActionStep: true` during config processing. It has no children, no links, no history. Designed as
a `visibleOn` target for actions from outer pipelines. The type guard `isPipelineActionConfig()` identifies
action configs before processing; after processing they are regular static pipeline configs.

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
