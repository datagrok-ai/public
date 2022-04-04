<!-- TITLE: Develop custom filters -->

# Filters

Developers can extend Datagrok with custom filters. This could be done by defining a class that extends
[DG.Filter](https://datagrok.ai/js-api/classes/dg.Filter) class. An
[example](https://github.com/datagrok-ai/public/blob/master/packages/Widgets/src/filters/radio-button-filter.ts) of such
class can be found in the [Widgets](https://github.com/datagrok-ai/public/tree/master/packages/Widgets) package.

The filter then must be registered in the `package.js` file:

```js
//name: Single Choice
//description: A filter that lets you select exactly one category
//tags: filter
//output: filter result
export function radioButtonFilter() {
  return new RadioButtonFilter();
}

```

The filter then can be invoked in the package with JS API as shown in
[custom filters example](https://dev.datagrok.ai/js/samples/ui/viewers/filters/custom-filters):

```js
let tv = grok.shell.addTableView(grok.data.demo.demog());
tv.filters({filters: [
  {type: 'Widgets:radioButtonFilter', columnName: 'race'},
]});

```

Or, alternatively, the filter can be added through UI:

![custom-filters](custom-filters.gif)

See also:

* [Filters](../../visualize/viewers/filters.md)
* [Filter class](https://datagrok.ai/js-api/classes/dg.Filter)
