# Document Schema: Datagrok Interactive Scientific Applications

## Common Part (reusable, identical for all applications)

```
datagrok-guide/
│
├── 1. Architecture and Philosophy
│   datagrok-interactive-app-guide.md
│   "Ports and adapters" pattern, pipelines, lifecycle,
│   subscriptions, teardown, testing.
│
├── 2. Datagrok API Reference
│   datagrok-api-reference.md
│   Catalog of inputs, viewers, buttons, tooltips, dialogs,
│   layouts, notifications, subscriptions. With links to documentation.
│
├── 3. Coding Standard
│   datagrok-coding-conventions.md
│   File structure, constants, types, classes, functions,
│   error handling, naming, formatting, testing.
│
├── 4. Specification Template
│   datagrok-app-specification-template.md
│   Empty template defining the specification structure
│   for a specific application.
│
├── 5. End-to-End Application Example
│   datagrok-example-app/
│   Completed specification + full code demonstrating
│   all architecture layers. Serves as a reference for generation.
│
└── 6. Implementation References (provided during core implementation phase)
    COMPUTATION-PATTERNS.md
    Working with raw data, null handling, single-pass aggregation,
    data locality, bool columns, module structure.

    ARRAY-OPERATIONS.md
    Typed array patterns: buffer reuse, out-parameter,
    scratch buffers, ring buffer, array pool, in-place transforms.

    WORKER-GUIDE.md
    Worker-utils infrastructure, WorkerColumn/WorkerDataFrame,
    DG↔Worker transforms, matrix layouts, lifecycle.

    PARALLEL-EXECUTION.md
    Fan-out/fan-in, distribution across workers, worker count.
```

## Application-Specific Part (unique for each application)

```
my-app-spec/
│
├── Application Specification (required)
│   my-app-specification.md
│   Completed specification template: core tasks, controls,
│   validation, reactivity, computation behavior,
│   rendering, layout, data lifecycle, UX.
│
├── Method Specifications (if the application has its own methods)
│   methods/
│   ├── method-A.md
│   │   Mathematical formulation, step-by-step algorithm,
│   │   inputs/outputs, constraints, edge cases,
│   │   literature references.
│   ├── method-B.md
│   └── ...
│
├── UI Component Specifications (if there are custom elements)
│   ui-components/
│   ├── component-X.md
│   │   Visual description (sketch/mockup), states,
│   │   events, styles (CSS), accessibility.
│   ├── component-Y.md
│   └── ...
│
└── External Library Documentation (if used)
    Not created, but referenced from the application specification.
    Links to API reference, README, guides.
```

## Relationship Between Parts

```
┌─────────────────────────────────────────────────────┐
│                  COMMON PART                         │
│                                                     │
│  ┌─────────────┐  ┌──────────────┐  ┌────────────┐ │
│  │Architecture │  │ API          │  │ Coding     │ │
│  │& Philosophy │  │ Reference    │  │ Standard   │ │
│  └─────────────┘  └──────────────┘  └────────────┘ │
│  ┌─────────────┐  ┌──────────────────────────────┐  │
│  │ Specification│  │  End-to-end example          │  │
│  │ Template    │  │  (specification + code)       │  │
│  └─────────────┘  └──────────────────────────────┘  │
│  ┌──────────────────────────────────────────────┐   │
│  │  Implementation References                   │   │
│  │  (computation, arrays, workers, parallelism) │   │
│  └──────────────────────────────────────────────┘   │
└──────────────────────────┬──────────────────────────┘
                           │
                           │ provided in context
                           │ together with
                           ▼
┌─────────────────────────────────────────────────────┐
│              APPLICATION-SPECIFIC PART               │
│                                                     │
│  ┌──────────────────────────────────────────────┐   │
│  │         Application Specification            │   │
│  │         (completed template)                  │   │
│  └──────┬──────────────┬───────────────┬────────┘   │
│         │              │               │            │
│         ▼              ▼               ▼            │
│  ┌────────────┐ ┌─────────────┐ ┌──────────────┐   │
│  │ Method     │ │ UI Component│ │ External     │   │
│  │ Specs      │ │ Specs       │ │ Library      │   │
│  │ (custom)   │ │ (custom)    │ │ Documentation│   │
│  └────────────┘ └─────────────┘ │ (links)      │   │
│                                 └──────────────┘   │
└─────────────────────────────────────────────────────┘
                           │
                           ▼
                  ┌─────────────────┐
                  │  Generated      │
                  │  Application    │
                  │  Code           │
                  └─────────────────┘
```

## Summary: What Is Needed to Create a New Application

| What | Source | Created from scratch? |
|---|---|---|
| Architecture and Philosophy | Common Part | No |
| Datagrok API Reference | Common Part | No |
| Coding Standard | Common Part | No |
| Specification Template | Common Part | No |
| End-to-End Example | Common Part | No |
| Implementation References | Common Part | No |
| Application Specification | Application-Specific Part | Yes, for each application |
| Method Specifications | Application-Specific Part | Yes, if the application has its own methods |
| UI Component Specifications | Application-Specific Part | Yes, if there are custom elements |
| External Library Documentation | Application-Specific Part (links) | No, already exists |
