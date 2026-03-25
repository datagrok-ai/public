# FlowCytometry Plugin for Datagrok — Agent Reference

This directory contains the complete reference documentation for building the `FlowCytometry` plugin for [Datagrok](https://datagrok.ai). It is structured for **agentic development**: each file is scoped to a specific concern so an agent can load only what it needs.

## Document Map

| File | Purpose | When to Load |
|---|---|---|
| `01-domain-overview.md` | What flow cytometry is, who uses it, the full analytical workflow | Understanding the problem domain |
| `02-file-formats.md` | FCS binary format spec (all versions), GatingML, workspace files, data concepts | Implementing parsers, importers, exporters |
| `03-example-files.md` | Canonical FCS test files with download URLs | Setting up test fixtures |
| `04-competitive-landscape.md` | All major competitors, pricing, strengths/weaknesses | Feature prioritization, positioning |
| `05-market-and-strategy.md` | Market size, growth drivers, Datagrok's commercial opportunity | Strategic decisions, monetization |
| `06-features.md` | Prioritized feature list by segment — MVP vs. advanced | Sprint planning, roadmap |
| `07-architecture.md` | Package structure, data model, Datagrok API integration points | All implementation tasks |
| `08-library-stack.md` | Every npm/Python/R library with URLs, rationale, usage pattern | Dependency selection, scaffolding |
| `09-workflows.md` | End-to-end user workflows: import → compensate → gate → analyze → export | UI/UX design, viewer design |
| `10-performance.md` | Strategies for large FCS files, Web Workers, WebGL, server-side compute | Performance-sensitive code |
| `11-roadmap.md` | Phased development plan with scope per phase | Sprint planning, prioritization |

## Plugin Identity

- **Package name:** `FlowCytometry`
- **Datagrok package type:** Plugin (TypeScript/Webpack)
- **Primary language:** TypeScript (browser) + Python (server-side scripting via Datagrok Conda)
- **Target users:** Pharma/biotech drug discovery immunologists, CRO scientists
- **Core value proposition:** Browser-native flow cytometry analysis integrated with Datagrok's compound registration, plate management, compliance, and scripting infrastructure

## Key Design Principles (for all agents)

1. **FCS files are binary** — always read as `ArrayBuffer`, parse with `DataView`. Never treat as text.
2. **Events are DataFrame rows** — each cell/event = one row; each channel/parameter = one column.
3. **Compensation before gating** — always apply compensation matrix before any gate evaluation.
4. **Gating is hierarchical** — gates form a tree; child gates operate only on parent-gated events.
5. **Logicle transform required** — raw and compensated data spans 5+ decades including negatives; log scale is wrong. Use logicle (biexponential) or arcsinh.
6. **Never block the main thread** — FCS parsing and dimensionality reduction go in Web Workers.
7. **Reuse Datagrok primitives** — use `DG.DataFrame`, `DG.JsViewer`, file import API, scripting API rather than building standalone.
8. **GatingML 2.0 for interop** — all gate definitions must be exportable as valid GatingML 2.0 XML.
