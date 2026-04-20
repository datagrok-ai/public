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
- Tree state (static/parallel/sequential pipeline nodes + FuncCall leaves)
- Data/validator/meta links between steps
- Consistency tracking and validation
- Serialization/deserialization of pipeline state

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

1. Rebuild compute-utils: `cd libraries/compute-utils && npm run build`
2. Rebuild Compute2: `cd packages/Compute2 && npm run build`
3. Publish: `cd packages/Compute2 && grok publish --release`
4. Run tests: see `packages/Compute2/CLAUDE.md` for full test instructions
