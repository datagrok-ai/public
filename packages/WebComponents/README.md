# WebComponents

WebComponents is a [package](https://datagrok.ai/help/develop/develop#packages) for the
[Datagrok](https://datagrok.ai) platform that provides runtime registration of
[Web Components](https://developer.mozilla.org/en-US/docs/Web/API/Web_components) as
custom HTML elements. It bridges the `@datagrok-libraries/webcomponents` library and the
platform, making reusable UI elements available throughout Datagrok.

## Registered Components

The package registers the following custom HTML elements on initialization:

| Element | Description |
|---------|-------------|
| `<dg-viewer>` | Viewer web component |
| `<dg-input-form>` | Input form web component |
| `<dg-button>` | Custom button (extends `<button>`) |
| `<dg-big-button>` | Larger button variant (extends `<button>`) |
| `<dg-icon-fa>` | Font Awesome icon |
| `<dg-icon-image>` | Image-based icon |
| `<dg-toggle-input>` | Toggle/checkbox input |
| `<dg-combo-popup>` | Combo box with popup dropdown |
| `<dg-markdown>` | Markdown renderer |
| `<dg-validation-icon>` | Validation status indicator |
| `<dock-spawn-ts>` | Docking panel manager (from `@datagrok-libraries/dock-spawn-dg`) |

All component implementations live in the
[`@datagrok-libraries/webcomponents`](../../libraries/webcomponents) library. This package
handles registration via `customElements.define()` so they are available at runtime.

## Package Structure

```
WebComponents/
  src/
    package.ts              # Entry point, registers all web components
    package-test.ts         # Test entry point
    styles.css              # Custom CSS for component styling
  detectors.js              # Semantic type detectors (placeholder)
  webpack.config.js
  package.json
```

## Build

```bash
npm install
npm run build              # grok api && grok check --soft && webpack
npm run lint               # ESLint check
npm run lint-fix           # ESLint auto-fix
```

## Development

Link local dependencies for development:

```bash
npm run link-all           # Links datagrok-api, webcomponents, and dock-spawn-dg
```

## See also

- [`@datagrok-libraries/webcomponents`](../../libraries/webcomponents) - Component implementations
- [`@datagrok-libraries/dock-spawn-dg`](../../libraries/dock-spawn-dg) - Docking manager
- [Datagrok JavaScript API](https://datagrok.ai/help/develop/packages/js-api)
- [Packages](https://datagrok.ai/help/develop/#packages)
