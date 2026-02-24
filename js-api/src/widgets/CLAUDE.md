# JS API Widgets - Input System

## Input Class Hierarchy

- `InputBase<T>` (`inputs-base.ts`) - wraps Dart-side inputs, provides events and property binding
- `JsInputBase<T>` (`inputs-base.ts`) - abstract base for pure-JavaScript inputs (no Dart UI)
- `MarkdownInput` (`markdown-input.ts`) - rich text editor using Quill.js

## Creating a Custom JS Input

Extend `JsInputBase<T>` and implement these 7 abstract members:

| Member | Purpose |
|--------|---------|
| `get inputType()` | Unique type name (e.g. `'Cron'`) |
| `get dataType()` | `DG.TYPE.*` constant for the underlying value |
| `getInput()` | Returns the root `HTMLElement` for the editor UI |
| `getValue()` | Returns the current typed value |
| `setValue(value)` | Sets value programmatically (fire `fireChanged()`) |
| `getStringValue()` | Returns the value as a string |
| `setStringValue(value)` | Sets the value from a string |

## Registration

Register in a package's `PackageFunctions` class:

```typescript
@grok.decorators.func({
  meta: {
    propertyType: 'string',   // matches DG.Property.propertyType
    semType: 'yourSemType',   // semantic type this editor handles
    role: 'valueEditor',      // tells the platform this is a value editor
  },
  outputs: [{type: 'object', name: 'result'}],
})
static yourInput(): DG.InputBase {
  return new YourInput();
}
```

The platform automatically uses this editor for properties with matching `semType`.

## Event Pattern

- `fireInput()` - call when the user interacts (typing, selecting). Triggers `onInput`.
- `fireChanged()` - call whenever the value changes (user or programmatic). Triggers `onChanged`.
- For user actions: fire both `fireInput()` then `fireChanged()`.
- For programmatic `setValue()`: fire only `fireChanged()`.

## Validation

Call `addValidator(fn)` in the constructor. The function receives the current value and returns `null` (valid) or an error message string.

## Reference Examples

- `FooInput` in `packages/Widgets/src/package.ts` - minimal two-field string input
- `MarkdownInput` in `js-api/src/widgets/markdown-input.ts` - rich text editor
- `CronInput` in `packages/PowerPack/src/widgets/cron-input.ts` - cron schedule editor with presets
