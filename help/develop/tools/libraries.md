<!-- TITLE: Library tour -->

# Library tour

This tour takes you through common libraries used to develop applications on top of the platform. Note that the list is
not comprehensive, and other third-party tools may be used along with the ones we mention.

## Datagrok toolkit

First of all, to work with your project, you need a utility to publish
[packages](../develop.md#packages) to the platform — `datagrok-tools`. See its
[documentation](https://github.com/datagrok-ai/public/tree/master/tools#datagrok-tools)
for installation instructions and usage examples. Upon completing package setup, you will see a set of default
dependencies listed in the `package.json`
file. One of them is `datagrok-api` that provides TypeScript API to the platform's core functionality. Whenever you need
details about a particular class or endpoint, consult the [API reference](https://datagrok.ai/js-api).

_Links:_

- [datagrok-tools](https://www.npmjs.com/package/datagrok-tools)
- [datagrok-api](https://www.npmjs.com/package/datagrok-api)

## Bundlers

Package templates come with a configuration for the `Webpack` bundler. In
`package.json`, `webpack` and `webpack-cli` are added as default development dependencies to build your package. When
working with various file extensions, you may have to include special modules
called [loaders](https://webpack.js.org/concepts/loaders). For better development experience, you can also
install `webpack-dev-server`.

_Links:_

- [webpack](https://www.npmjs.com/package/webpack)
- [webpack-cli](https://www.npmjs.com/package/webpack-cli)
- [webpack-dev-server](https://www.npmjs.com/package/webpack-dev-server)

## TypeScript

As we recommend [TypeScript](https://www.typescriptlang.org/) as a language for package development, there's an option
to [create a package](../develop.md#getting-started)
with a `--ts` flag. Among other things, it adds two new dependencies: `typescript` (provides the language support)
and `ts-loader` (a file loader for `webpack`). Likely, these are not the only libraries you will use when writing in
TypeScript. For example, if one of the libraries you want to work with has not been typed yet, check out the
[Definitely Typed](https://github.com/DefinitelyTyped/DefinitelyTyped) resource, which provides type definitions for
popular packages.

_Links:_

- [TypeScript](https://www.npmjs.com/package/typescript)
- [ts-loader](https://www.npmjs.com/package/ts-loader)

## Visualization

The platform comes with a diverse set of visualizations (see the
[Viewers](../../visualize/viewers.md) article). Moreover, it is possible to
[create a custom viewer](../how-to/develop-custom-viewer.md) using our API. For this task, you can use such libraries
as `d3`, `three.js`, or `echarts`. Datagrok's public repository contains packages with examples:
[Viewers](https://github.com/datagrok-ai/public/tree/master/packages/Viewers)
and [Charts](https://github.com/datagrok-ai/public/tree/master/packages/Charts). For scientific applications, you may
find useful such projects as
[NGL](https://www.npmjs.com/package/ngl), [RDKit](https://www.npmjs.com/package/@rdkit/rdkit),
and [OpenChemLib](https://www.npmjs.com/package/openchemlib), but first look at what solutions already exist to
integrate with them (see the [Cheminformatics](../domains/chem/cheminformatics.md)
page).

_Links:_

- [d3](https://www.npmjs.com/package/d3)
- [three](https://www.npmjs.com/package/three)
- [echarts](https://www.npmjs.com/package/echarts)
