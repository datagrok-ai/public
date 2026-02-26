# Compute2

Compute2 is a [package](https://datagrok.ai/help/develop/develop#packages) for the
[Datagrok](https://datagrok.ai) platform that provides the modeling UI and workflow
execution framework. It powers the [Compute](https://datagrok.ai/help/compute/)
capabilities of the platform, including rich function editors, multi-step workflow
pipelines, and the Model Hub.

## Features

### Rich Function View (RFV)

The primary editor for executing individual Datagrok functions with a full-featured UI:

- Input form with validation and consistency checking
- Multi-tab output visualization for scalars and dataframes
- Run history browser with comparison support
- Sensitivity analysis
- Fitting and optimization (including formula-based fitting)
- Help panel with function documentation
- Export to Excel with custom formatting

### Tree Wizard (Workflows)

A visual editor for building and executing sequential and parallel workflow pipelines:

- Tree navigation with step status indicators
- Support for static, sequential, and parallel pipelines
- Drag-and-drop step reordering
- Dynamic step addition and removal
- Input/output linking between steps with validation
- Metadata handlers for dynamic configuration
- Subtree execution and consistency management
- Full workflow export to ZIP with Excel sheets
- Buttons-based navigation mode

### Model Hub

An application for browsing and launching computational models, accessible from the
Compute menu in the platform ribbon.

### Custom Function Views

Extensible base class (`CustomFunctionView`) for building custom UIs for function
execution beyond what the default RFV provides.

## Technology

Built with Vue 3 (Composition API), RxJS for reactive state management, and Tailwind CSS
for styling. Uses the `@datagrok-libraries/compute-utils` library for workflow computation
and history management.

## Demo Scripts

The package includes server-side demo scripts in `scripts/`:

- **ObjectCooling2** - Newton's law of cooling simulation with sensitivity analysis and fitting metadata
- **LongScript** - Long-running script simulation with multiple outputs
- **SimpleInputs** - Basic input type demonstration

## Package Structure

```
Compute2/
  src/
    package.ts              # Entry point, function and editor registration
    apps/                   # Standalone Vue applications (RFV, TreeWizard, History)
    components/             # Reusable Vue components
      RFV/                  # Rich Function View components
      TreeWizard/           # Workflow tree navigation
      PipelineView/         # Pipeline flow visualization
      History/              # Run history browser
      Inspector/            # Developer debug tool
      Logger/               # Real-time log viewer
    composables/            # Vue 3 composition utilities
    test/                   # Test utilities
  scripts/                  # Server-side demo scripts (JS)
  files/                    # Static assets (icons)
  webpack.config.js
  package.json
```

## Build

```bash
npm install
npm run build              # grok api && grok check --soft && webpack
npm run build-all          # Build js-api, libraries, and this package in order
npm run lint               # ESLint check
npm run lint-fix           # ESLint auto-fix
```

## Development

Link local dependencies for development:

```bash
npm run link-all           # Links datagrok-api and @datagrok-libraries/*
```

## See also

- [Compute documentation](https://datagrok.ai/help/compute/)
- [Datagrok JavaScript API](https://datagrok.ai/help/develop/packages/js-api)
- [Packages](https://datagrok.ai/help/develop/#packages)
