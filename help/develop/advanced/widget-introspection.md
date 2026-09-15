---
title: "Widget introspection"
description: How AI agents, command palettes and UI tests discover what is on the screen — widget descriptors, runtime status, AI briefings, functions, and the property catalogs of platform entities.
keywords:
  - WidgetDescriptor
  - getWidgetStatus
  - aiDescription
  - getFunctions
  - grok.meta
  - AI agents
  - UI testing
---

Every widget on the screen — viewers, views, inputs, dialogs — answers the same small set of
questions without widget-specific code. The AI assistant, command palettes and UI tests are built
on these answers, and a plugin gets them for free by extending `DG.Widget`.

| Question                                   | API                                                                                                     |
|--------------------------------------------|---------------------------------------------------------------------------------------------------------|
| Which widgets exist, without creating one? | `DG.WidgetDescriptor.getDescriptors()`, `getByName(name)`: name, synonyms, description, properties, events, icon |
| What is on the screen now?                 | `DG.Widget.getAll()`; each widget's `descriptor`                                                        |
| What state is it in?                       | `widget.getWidgetStatus()`: `parts`, `hitAreas`, `shortcuts`, `events`, `description`, `error`, `inputs` |
| What can be done with it?                  | `widget.getFunctions()` (real, typed `DG.Func`s) and `widget.props` (the property bag)                 |
| How should an assistant approach it?       | `widget.aiDescription`, a short briefing shown to the assistant as workspace context                     |
| What are the fields of a platform type?    | `grok.meta.propertiesOf(type)`, `grok.meta.coreLocationOf(type)`                                        |

## Descriptors

A descriptor is the static half: what a widget type is called, what it is for, which properties
and events it has. It is available before any instance exists, which is what a command palette
or a "which viewer fits this question" prompt needs.

```js
for (const d of DG.WidgetDescriptor.getDescriptors())
  console.log(d.name, d.synonyms, d.properties.map((p) => p.name), d.events.map((e) => e.eventName));

const sp = DG.WidgetDescriptor.getByName(DG.VIEWER.SCATTER_PLOT);
```

Live sample: [descriptors](https://public.datagrok.ai/js/samples/ui/viewers/advanced/descriptors).

## Runtime status

`getWidgetStatus()` is the dynamic half: a snapshot a test or an agent can act on.

* `parts`: named DOM elements (`overlay`, `xAxis`, a form's fields) to locate and interact with.
* `hitAreas`: named rectangles to click at, in the widget's own coordinates.
* `shortcuts`: key combination to a human-readable description.
* `events`: what the widget fires; subscribe with `widget.onEvent(name)`.
* `description`: free-form summary of the current state (what is displayed, how many rows).
* `error`: the validation message, or null when the widget is in a valid state.
* `inputs`: for forms, one record per field with its value, choices, validity and error.

```js
const sp = view.scatterPlot({x: 'height', y: 'weight'});
const s = sp.getWidgetStatus();
for (const el of Object.values(s.parts))
  el.style.outline = '2px solid red';
sp.onEvent('d4-scatterplot-point-click').subscribe((e) => grok.shell.info(`row ${e.args.rowId}`));
```

Dart viewers override the status with viewer-specific parts and hit areas; a JS widget overrides
`getWidgetStatus()` to expose its own.

Live samples: [widget status](https://public.datagrok.ai/js/samples/ui/interactivity/widget-status),
[the machine surface of a form](https://public.datagrok.ai/js/samples/dapi/domains/facade-machine-surface).

## Functions and properties

Actions are registered functions, not ad-hoc methods: `getFunctions()` returns the `DG.Func`s
applicable to the widget, with typed parameters and descriptions, and the assistant calls them the
same way the context menu does. Named values are properties: writing `widget.props.priority =
'critical'` takes the same path as typing into the field, including validation, so the next
`getWidgetStatus().error` reflects it.

```js
const w = DG.Widget.getAll().find((w) => w.type === 'DomainForm');
w.props['priority'] = 'critical';
const save = w.getFunctions().find((f) => f.name === 'Save');
await save.apply({widget: w});
```

## AI briefing

`aiDescription` is a few sentences for the assistant: what the widget is, what its functions do,
and where to start. Platform widgets ship a default (the console, the sketcher, the browse panel,
the permissions browser); u2 components seed it from their registry usage or description. Set it on
your own widgets, views and dialogs where the default would leave the assistant guessing:

```js
const view = DG.View.create();
view.name = 'Plate reader';
view.aiDescription = 'Shows the last uploaded plate. Call "Normalize" before any statistics; ' +
  '"Export" writes the normalized plate as CSV.';
```

## Entity properties

`grok.meta` answers the same question for data rather than widgets: which properties a platform
entity type (`'User'`, `'Project'`) or a domain table (`'Core.users'`, `'plates.plate'`) exposes.
The names it returns are the names the platform's grids, filter panels and server-side facets use,
so discovery and enforcement do not drift apart.

```js
const props = await grok.meta.propertiesOf('User', {filterable: true});
const users = await grok.dapi.users.filter(`${props[0].name} = "admin"`).list();
const where = await grok.meta.coreLocationOf('User');   // {schema: 'Core', table: 'users'} or null
```

Live sample: [entity properties](https://public.datagrok.ai/js/samples/dapi/entity-properties).

See also:

* [JS API reference](https://datagrok.ai/api/js): `WidgetDescriptor`, `Widget`, `Viewer`, `Meta`
* [Manipulating viewers](../how-to/viewers/manipulate-viewers.md)
* [Datagrok UI](ui.md)
