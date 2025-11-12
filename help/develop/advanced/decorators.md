---
title: "Decorators"
---

# Decorators

Package functions are typically registered in the main package file `package.ts`. Each function has a special [parameter annotation](../../datagrok/concepts/functions/func-params-annotation.md) depending on its [role](../function-roles.md).
You can use [decorators](https://www.typescriptlang.org/docs/handbook/decorators.html) to register annotated functions based on decorated class/object Decorators add the ability to create strong type annotation instead of writing it by hand.

You can use a decorator `@grok.decorators.<name>`. This is equivalent to adding a function to `package.ts`. When you run the `build` script for your package, the webpack plugin called `FuncGeneratorPlugin` will add a special `package.g.ts` file to your project. Note that it is not on the ignore list, so you should commit this file to the repository.

## Class decorators

There are a few [Datagrok JS API](https://datagrok.ai/api/js) that can be decorated to register them as annotated functions. You are also able to decorate their subclasses. Here is the list of classes for which decorators are designed for: 

* [custom viewers](../how-to/viewers/develop-custom-viewer.md)
* [custom filters](../how-to/viewers/custom-filters.md)
* [custom cell renderers](../how-to/grid/custom-cell-renderers.md)

### Examples

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

## Method decorators

Additionally, it is possible to decorate methods. There is no ability to decorate functions, which is why you should decorate static methods in any class you need. The best class for that is `PackageFunctions` in the `package.ts` file.

This functionality enables the automatic generation of annotated functions from decorated class methods. Base annotation gets created from the method's signature, which is an option in the decorator's input options. You can see [list](https://github.com/datagrok-ai/public/blob/master/js-api/src/decorators/functions.ts) of decorators and their options in the js-api package.

There is`@grok.decorators.param` decorator that allows setting input options. You can easily set the name, type, and options in the decorator's input object. Most of the types are already [registered](https://github.com/datagrok-ai/public/blob/7ce4edaca2a937b934120a8894ea927c1eca8c7f/js-api/src/const.ts#L82) in `DG.TYPE` enum and you can easily use it to set params type.

### Example 

<details>
<summary> Func Decorators </summary>
<div>

```ts
export class PackageFunctions { 

  @grok.decorators.init({})
  static async init() {}
  
  @grok.decorators.func({ name: 'solveODE'})
  static solve(problem: ODEs): DG.DataFrame {
    return solveDefault(problem);
  }
}
```

</div>
</details>



<details>
<summary> Params Decorator</summary>
<div>

```ts
export class PackageFunctions { 
  @grok.decorators.func({
    name: 'Ball flight',
    description: 'Ball flight simulation',
    meta: {
      icon: 'files/icons/ball.png'
    },
    outputs: [
      {
        name: 'maxDist', 
        type: 'double',
        options: {caption: 'Max distance'}
      }
    ]
  })
  static ballFlight(
    @grok.decorators.param({options:{initialValue: '0.01', category: 'Ball', caption: 'Diameter', units: 'm', min: '0.01', max: '0.3'}}) dB: number, 
    @grok.decorators.param({options:{initialValue: '200', category: 'Ball', caption: 'Density', description: 'Material density', units: 'kg/m^3', min: '200', max: '1200'}}) roB: number, 
    @grok.decorators.param({type: DG.TYPE.INT, options:{initialValue: '50', category: 'Throw parameters', caption: 'Velocity', min: '40', max: '60', units: 'm/sec'}}) v: number, 
    @grok.decorators.param({type: DG.TYPE.FLOAT, options:{initialValue: '45', category: 'Throw parameters', caption: 'Angle', min: '20', max: '70', units: 'deg'}}) a: number) {
  }
}
```

</div>
</details>


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
     "datagrok-tools": "^4.14.1"
   }
   ```

1. After the build ensure that there is the next line in the `package.ts` exists(if it doesnt add it):

   ```ts
   export * from './package.g';
   ```

See also:

* [Function roles](../function-roles.md)
* [Parameter annotation](../../datagrok/concepts/functions/func-params-annotation.md)
* [JavaScript development](../develop.md)
* [TypeScript class decorators](https://www.typescriptlang.org/docs/handbook/decorators.html#class-decorators)
