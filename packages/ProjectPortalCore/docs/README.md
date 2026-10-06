# Project Portal Core: design docs

The package [README](../README.md) is the usage guide. These pages hold the design behind it.

| Page | What it answers |
|---|---|
| [GOAL.md](GOAL.md) | What the portal is for, in 15 points, and who delivers each: the core, a customer package or the platform |
| [multi-portfolio-design.md](multi-portfolio-design.md) | Why the portfolio is the unit of configuration, and what real registries raise when several portfolios share one instance |
| [architecture.md](architecture.md) | The reference architecture: layers, the core/customer package split, schema design rules, the workspace builder contract, and where the build differs |
| [features.md](features.md) | How each feature behaves, in detail |
| [schema.md](schema.md) | The `portal` schema: table groups, who may write each table, choices made while building, migrations, version history |
| [build-log.md](build-log.md) | Each feature's schema change and the decisions taken on the way that are open to revisit |
| [platform-findings.md](platform-findings.md) | Platform issues found while building, their workarounds, and what the server enforces on writes |
| [testing.md](testing.md) | What the package tests and BDD scenarios cover, and the known failures |

Ids such as F-12, FR-50, PORT-35 or PLAT-10 are the product's feature, requirement, story and platform-ask ids,
kept so decisions can be traced; GROK-nnnnn are platform tickets. Every example uses the synthetic portfolios the
package ships (EX-SMD, EX-BIOD, EX-BIOC); customer configurations live in their own packages, outside this
repository.
