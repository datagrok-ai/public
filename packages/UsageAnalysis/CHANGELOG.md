# Usage Analysis changelog

## v.next

* GROK-6166: Query viewers: Surface query-load failures in the UI — missing/failing queries now clear the loader spinner and show an error balloon instead of silently swallowing the rejection
* Stress tests: Replaced the old stress dashboards/queries with `StressTestsSummary` (all builds — median/avg/min/max/p95 of ms by threads, pass rate as the mean per-test pass fraction), `StressTestsRaw` (per-build raw runs incl. error text, latest default, boxplot of ms), and `StressTestsFailures` (failing test×thread combos with a sample error); rewired the CI/CD tests project. Backed by the new server `stress_tests` table populated from the Jenkins stress job.
* GROK-14456: Usage Analysis: Log tab improvements (added parameter details to the context panel and stack traces for errors)
* GROK-12108: Usage Analysis: Errors tab
* GROK-19820: Jira swagger: Removed three unused endpoints, keeping only `Jira Create Issue`
* Metrics: Added Admin Metrics dashboard tab
* TestTrack: Fixed node-expansion errors after the `Test Track`→`TestTrack` folder rename — the loader now skips `-run.md` run reports and `-spec.ts` files (filtered to `.md`) instead of parsing them as test cases, and guards against a null category node
* GROK-16262: TestTrack: Removed the `Browse > Browse tree states` test case

## 2.5.1 (2026-05-21)

* TestTrack: Wired Playwright suite (`files/TestTrack`) into `grok test` via the new `playwrightTests` opt-in; specs now authenticate via the same dev-key→token flow as the Puppeteer pass and run in CI alongside existing tests
* Test Dashboard: Moved "Watch a jira ticket" storage from the `ua_tickets` Postgres schema to a JSON file in `System:AppData/UsageAnalysis/manual-tickets.json`; removes the package-owned schema (and its `UA_tickets` connection / `ManualTicketFetch` / `ManualTicketCreation` queries) so UA installs cleanly on customer envs without CREATE SCHEMA/ROLE privileges
* Test Dashboard: Wired the JIRA-ticket verdict grids back into `onFrameAttached` so watched/auto-detected tickets now render as priority-bucketed grids in the dashboard accordion (previously the display path was defined but never called)
* Test Dashboard: Removed dead code (~115 lines): unused `verdictsOnTestTrack`/`testTrackDowngrades`, commented-out `unaddressedTests` block, stray `debugger` statement, and several unused locals
* Click Events widget: Added date range filter, section labels, highlight zones tooltip

## 2.5.0 (2026-03-20)

### Features

* GROK-18787: Click Tracking:
  * Added Clicks tab
  * Added tab control for Click Analysis and Specific Clicks
  * Added elements inspect
* [#3213](https://github.com/datagrok-ai/public/issues/3213): Introduced Projects tab
* Improved query performance, load time, and projects tab
* Reports: Added descending ordering by creation date in widget
* Reports: Truncated long text in widget
* Improved ticket loading with bulk method and background processing
* Decreased bundle size by removing unnecessary test library usages
* [#3604](https://github.com/datagrok-ai/public/issues/3604): Migrated tags to roles

### Bug Fixes

* GROK-19855: Reporting app now loads reports when opened by URL with report number
* GROK-19523: AppTreeBrowser decorator now finds the app properly
* GROK-19504: Invalid datagrok-api imports corrected

## 2.4.2 (2025-09-26)

* Test Track: Added advanced test case for visual query
* Reports: Fixed not visible add rule button
* TestTrack: Removed manual test that have autotest
* GROK-18772: ManualTests Dashboard: Unable to parse string "Test Track:test" when trying to "Run" test
* Test Track: updated and refined test cases for Projects (added Spaces)
* Tests: Added warning about missing packages tests
* Update Autocomplete.md
* Update Formula Refreshing.md
* GROK-17332: Removed unused queries
* Table name change
* Test Track: Line chart: Added context menu for test cases
* GROK-18366 Usage analysis: Log: an error occurs when opening the tab
* UsageAnalysis: Api: naming fixes

## 2.3.2 (2025-03-24)

* Tests Dashboards: improvements for gathering of benchmarks

## 2.3.1 (2025-03-13)

* Dependency: datagarok-api >= 1.24.0*

### Features

* Tests Dashboards: improvements and bug fixes
* Stress tests dashboards revisited 

## 2.3.0 (2025-02-18)

* Dependency: datagarok-api >= 1.24.0*

### Features

* Packages: Possibility to filter packages by category
* Functions: Possibility to filter functions by tags
* Improvements in usage queries and bug fixes
* Tests Dashboard

## 2.2.1 (2024-12-03)

* Dependency: datagarok-api >= 1.22.0*

### Features:

* Added detector for tests tables
* Updated test track business logic
* Added test cases to Test Track app 

## 2.2.0 (2024-10-24)

* Dependency: datagarok-api >= 1.22.0*

### Features:

* Reporting: Better subscriptions handling
* Service logs app: Application that gives possibility to see logs of different services within Datagrok infrastructure.

## 2.0.2 (2024-07-23)

* Dependency: datagarok-api >= 1.20.0*

### Features

* Reporting: Bugs fixes and improvements in UX
* Test Track:
  * Collaborative testing synchronization
  * Ability to set severity level of errors
  * Ability to associate test executions with Error reports
  * Updates Test Track according new tests registration

## 2.0.1 (2024-06-12)

### Features

* Reporting app

## 2.0.0 (2024-03-16)

* Dependency: datagarok-api >= 1.18.0*

### Features

* Test Track released
* GROK-14344 UA: Functions execution time
* UA: update routing, overall improvements
* UA: Test results for selected date in Test Track
* UA: Compare test results in Test Track
* GROK-14402 UA: Test Track

## 1.0.4 (2023-11-27)

*Dependency: datagarok-api >= 1.17.3*

### Features

* GROK-13271 Usage Analysis: Logs view (#2363)
* UA: Refactor Usage tab for Packages & fixes

### Bug Fixes

* UA: Tests fixes
* Enabled cards caching

## 1.0.3 (2023-11-09)

*Dependency: datagarok-api >= 1.16.1*

### Features

* UA: Packages installation time
* Overall improvements

## 1.0.2 (2023-10-31)

*Dependency: datagarok-api >= 1.16.1*

### Features

* GROK-13425 UA: Ability to call query without loading from server
* GROK-12578 UA: Tests
* GROK-13721 Test Track (WIP)
* UA:
    * Speedup
    * CSS tweaks
    * Better choiceInput
    * Client side cache
* GROK-13947 UA | TM: Usage tab for packages

### Bug Fixes

* UA: Fixes GROK-13498 & GROK-13500
* UA: Queries & widget fixes

## 1.0.1 (2023-07-03)

*Dependency: datagarok-api >= 1.13.8*

### Features

* Speedup, cleanup and styles fixes
* GROK-13342 Resurrect test tracking system

### Bug Fixes

* Query deploy bugfix

## 1.0.0 (2023-05-10)

*Dependency: datagarok-api >= 1.13.8*

### Features

* GROK-12111 UsageAnalysis: Packages tab
* GROK-12745 UsageAnalysis: Functions tab
* GROK-12796 UsageAnalysis: Widget
* GROK-12110: Usage Analysis: Overview tab
* GROK-12109 UsageAnalysis: Events tab
* UsageAnalysis: Cards
* Overall fixes & improvements

### Bug Fixes

* GROK-13127 UsageAnalysis: error in widget