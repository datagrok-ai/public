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

## After Making Changes

1. Rebuild compute-utils: `cd libraries/compute-utils && npm run build`
2. Rebuild Compute2: `cd packages/Compute2 && npm run build`
3. Publish: `cd packages/Compute2 && grok publish --release`
4. Run tests: see `packages/Compute2/CLAUDE.md` for full test instructions
