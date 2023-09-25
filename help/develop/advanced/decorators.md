---
title: "Class decorators"
---

# Class decorators

Package functions are typically registered in the main package file `package.ts`. Each function has a special
[parameter annotation](../../datagrok/concepts/functions/func-params-annotation.md) depending on its [role](../function-roles.md).
You can use [class decorators](https://www.typescriptlang.org/docs/handbook/decorators.html#class-decorators) to register
such package functions as:

* [custom viewers](../how-to/develop-custom-viewer.md)
* [custom filters](../how-to/custom-filters.md)
* [custom cell renderers](../how-to/custom-cell-renderers.md)

If a function uses a subclass that extends classes from [Datagrok JS API](https://datagrok.ai/js-api), you can use a
decorator `@grok.decorators.<name>`. This is equivalent to adding a function to `package.ts`. There is no need
to add anything other than the class itself. When you run the `build` script for your package, the webpack plugin called
`FuncGeneratorPlugin` will add a special `package.g.ts` file to your project. Note that it is not on the ignore list, so
you should commit this file to the repository.

## Enabling decorators in a package

1. In `tsconfig.json`, enable these experimental options:

   ```json
   "experimentalDecorators": true,   /* Enables experimental support for ES7 decorators. */
   "emitDecoratorMetadata": true,    /* Enables experimental support for emitting type metadata for decorators. */
   ```

1. In `webpack.config.js`, enable `FuncGeneratorPlugin`:

   ```js
   const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');

   module.exports = {
     /** ... */
     plugins: [
       new FuncGeneratorPlugin({outputPath: './src/package.g.ts'}),
     ],
   };
   ```

1. In `package.json`, add `datagrok-tools` to dev dependencies; upgrade its version, if necessary:

   ```json
   "devDependencies": {
     "datagrok-tools": "^4.12.x"
   }
   ```

## Examples

<details>
<summary> Custom viewers </summary>
<div>

```ts
@grok.decorators.viewer({
  name: 'Test Viewer',
  description: 'Creates a Test Viewer instance',
  icon: 'images/icon.png',
  toolbox: true,
  viewerPath: 'Tests | Show results',
})
export class TestViewer extends DG.JsViewer {
}
```

</div>
</details>

<details>
<summary> Custom filters </summary>
<div>

```ts
@grok.decorators.filter({
  name: 'Radio Button Filter',
  description: 'Single option filter',
  semType: 'Country',
})
export class RadioButtonFilter extends DG.Filter {
}
```

</div>
</details>

<details>
<summary> Custom cell renderers </summary>
<div>

```ts
@grok.decorators.cellRenderer({
  name: 'Fasta Sequence Cell Renderer',
  description: 'Macromolecule renderer',
  cellType: 'sequence',
  columnTags: 'quality=Macromolecule, units=fasta',
})
export class MacromoleculeSequenceCellRenderer extends DG.GridCellRenderer {
}
```

</div>
</details>

See also:

* [Function roles](../function-roles.md)
* [Parameter annotation](../../datagrok/concepts/functions/func-params-annotation.md)
* [JavaScript development](../develop.md)
* [TypeScript class decorators](https://www.typescriptlang.org/docs/handbook/decorators.html#class-decorators)
