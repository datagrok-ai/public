---
title: "Library tour"
description: A tour of common libraries used to develop Datagrok applications, including bundlers, TypeScript, and visualization tools.
keywords:
  - datagrok-tools
  - datagrok-api
  - rspack
  - TypeScript
  - d3
  - three.js
  - echarts
  - custom viewer libraries
---

This tour takes you through common libraries used to develop applications on top of the platform. Note that the list is
not comprehensive, and other third-party tools may be used along with the ones we mention.

## Datagrok toolkit

First of all, to work with your project, you need a utility to publish
[packages](../../develop.md#packages) to the platform — `datagrok-tools`. See its
[documentation](https://github.com/datagrok-ai/public/tree/master/tools#datagrok-tools)
for installation instructions and usage examples. Upon completing package setup, you will see a set of default
dependencies listed in the `package.json`
file. One of them is `datagrok-api` that provides TypeScript API to the platform's core functionality. Whenever you need
details about a particular class or endpoint, consult the [API reference](https://datagrok.ai/api/js).

_Links:_

- [datagrok-tools](https://www.npmjs.com/package/datagrok-tools)
- [datagrok-api](https://www.npmjs.com/package/datagrok-api)

## Bundlers

Packages are bundled by [rspack](https://rspack.rs) with swc, through the single `@datagrok/build-config`
development dependency: the configuration is shared by every Datagrok package, and a package adds an
`rspack.config.js` only for what differs (extra externals, WebAssembly, JSX). CSS, images and `.wasm` files
are handled out of the box; see the
[@datagrok/build-config README](https://github.com/datagrok-ai/public/blob/master/build-config/README.md).

_Links:_

- [@datagrok/build-config](https://www.npmjs.com/package/@datagrok/build-config)
- [rspack](https://www.npmjs.com/package/@rspack/core)

## TypeScript

As we recommend [TypeScript](https://www.typescriptlang.org/) as a language for package development, there's an option
to [create a package](../../onboarding/getting-started.md)
with a `--ts` flag. TypeScript itself comes with `@datagrok/build-config`: swc transpiles the sources when
bundling, and `npm run typecheck` runs the TypeScript compiler (`tsc --noEmit`) separately. If one of the
libraries you want to work with has not been typed yet, check out the
[Definitely Typed](https://github.com/DefinitelyTyped/DefinitelyTyped) resource, which provides type definitions for
popular packages.

_Links:_

- [TypeScript](https://www.npmjs.com/package/typescript)

## Visualization

The platform comes with a diverse set of visualizations (see the
[Viewers](../../../visualize/viewers/viewers.md) article). Moreover, it is possible to
[create a custom viewer](../../how-to/viewers/develop-custom-viewer.md) using our API. For this task, you can use such libraries
as `d3`, `three.js`, or `echarts`. Datagrok’s public repository offers an example implementation in the [Charts](https://github.com/datagrok-ai/public/tree/master/packages/Charts) package. For scientific applications, you may
find useful such projects as
[NGL](https://www.npmjs.com/package/ngl), [RDKit](https://www.npmjs.com/package/@rdkit/rdkit),
and [OpenChemLib](https://www.npmjs.com/package/openchemlib), but first look at what solutions already exist to
integrate with them (see the [Cheminformatics](../../../datagrok/solutions/domains/chem/chem.md)
page).

_Links:_

- [d3](https://www.npmjs.com/package/d3)
- [three](https://www.npmjs.com/package/three)
- [echarts](https://www.npmjs.com/package/echarts)
