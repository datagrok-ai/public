/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/sankey-chord/sankey-chord-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.sankey, charts.viewer.chord, public#4420, public#4421]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterBetween, filterPasses, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, propertyShouldBe, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Sankey and Chord edge cases — Source equal to Target, too many categories, a removed column", () => {
  const session = feature(test, "features/viewers/sankey-chord/sankey-chord-edge-cases.feature", import.meta.url);
  test("A Sankey whose Source and Target are the same column says the graph has cycles", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "When user adds a sankey viewer", () => addViewer(page, "sankey"));
    await session.step(20, "Then \"Source\" property of sankey viewer should be \"SEX\"", () => propertyShouldBe(page, "Source", el("sankey viewer"), "SEX"));
    await session.step(21, "When user sets \"Target\" property of sankey viewer to \"SEX\"", () => setProperty(page, "Target", el("sankey viewer"), "SEX"));
    await session.step(22, "Then \"The graph contains cycles. Please remove circular dependencies.\" text in sankey viewer should be visible", () => shouldBe(page, el("\"The graph contains cycles. Please remove circular dependencies.\" text in sankey viewer"), "visible"));
    await session.step(23, "And no errors should have been logged", () => noErrors(page));
    await session.step(24, "When user sets \"Target\" property of sankey viewer to \"RACE\"", () => setProperty(page, "Target", el("sankey viewer"), "RACE"));
    await session.step(25, "Then \"Caucasian\" text in sankey viewer should be visible", () => shouldBe(page, el("\"Caucasian\" text in sankey viewer"), "visible"));
    await session.step(26, "And \"The graph contains cycles. Please remove circular dependencies.\" text in sankey viewer should be absent", () => shouldBe(page, el("\"The graph contains cycles. Please remove circular dependencies.\" text in sankey viewer"), "absent"));
    await session.step(27, "And no errors should have been logged", () => noErrors(page));
  });
  test("A Chord over more than 1500 category pairs says there are too many categories", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(30, "When user adds a chord viewer", () => addViewer(page, "chord"));
    await session.step(31, "Then \"From\" property of chord viewer should be \"SEX\"", () => propertyShouldBe(page, "From", el("chord viewer"), "SEX"));
    await session.step(32, "When user sets \"To\" property of chord viewer to \"USUBJID\"", () => setProperty(page, "To", el("chord viewer"), "USUBJID"));
    await session.step(33, "Then \"Too many categories to render.\" text in chord viewer should be visible", () => shouldBe(page, el("\"Too many categories to render.\" text in chord viewer"), "visible"));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
    await session.step(35, "When user sets \"To\" property of chord viewer to \"RACE\"", () => setProperty(page, "To", el("chord viewer"), "RACE"));
    await session.step(36, "Then \"Asian\" text in chord viewer should be visible", () => shouldBe(page, el("\"Asian\" text in chord viewer"), "visible"));
    await session.step(37, "And \"Too many categories to render.\" text in chord viewer should be absent", () => shouldBe(page, el("\"Too many categories to render.\" text in chord viewer"), "absent"));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("The Sankey's Source column removed from the table raises no error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord", "@known-failure:50", "@realizes:public#4420"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(44, "When user adds a sankey viewer", () => addViewer(page, "sankey"));
      await session.step(45, "Then \"Source\" property of sankey viewer should be \"SEX\"", () => propertyShouldBe(page, "Source", el("sankey viewer"), "SEX"));
      await session.step(46, "When user removes \"SEX\" column", () => removeColumn(page, "SEX"));
      await session.step(47, "Then the table should not have a column \"SEX\"", () => hasNoColumn(page, "SEX"));
      await session.step(48, "When user filters rows where \"AGE\" is between 30 and 50", () => filterBetween(page, "AGE", 30, 50));
      await session.step(49, "Then 494 rows should pass the filter", () => filterPasses(page, 494));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    }, 50);
  });
  test("The Chord's From column removed from the table raises no error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord", "@known-failure:62", "@realizes:public#4421"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await knownFailure(async () => {
      await session.step(56, "When user adds a chord viewer", () => addViewer(page, "chord"));
      await session.step(57, "Then \"From\" property of chord viewer should be \"SEX\"", () => propertyShouldBe(page, "From", el("chord viewer"), "SEX"));
      await session.step(58, "When user removes \"SEX\" column", () => removeColumn(page, "SEX"));
      await session.step(59, "Then the table should not have a column \"SEX\"", () => hasNoColumn(page, "SEX"));
      await session.step(60, "When user filters rows where \"AGE\" is between 30 and 50", () => filterBetween(page, "AGE", 30, 50));
      await session.step(61, "Then 494 rows should pass the filter", () => filterPasses(page, 494));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    }, 62);
  });
});
