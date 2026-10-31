/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/radar/radar-current-row.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.radar, public#1472]
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
import {check, clickOn, isExpanded, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaColor, areaNotColor, clickArea, noErrors, propertyShouldBe, repainted, setProperty, takeSnapshot, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Radar: Show Current Row, and Pick Up / Apply between two Radars", () => {
  const session = feature(test, "features/viewers/radar/radar-current-row.feature", import.meta.url);
  test("Radar: Show Current Row, and Pick Up / Apply between two Radars", {tag: ["@journey", "@viewers", "@realizes:charts.viewer.radar", "@realizes:public#1472"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(19, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(20, "When user clicks on first \"Radar\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Radar\" card in \"Add Viewer\" dialog")));
    await session.step(21, "Then radar viewer should be visible", () => shouldBe(page, el("radar viewer"), "visible"));
    await session.step(22, "And \"Values\" property of radar viewer should be \"AGE, HEIGHT, WEIGHT\"", () => propertyShouldBe(page, "Values", el("radar viewer"), "AGE, HEIGHT, WEIGHT"));
    await run.scenario("Show Current Row draws the current row's line in green, and off it is gone", async () => {
      await session.step(25, "When user sets \"Show Mouse Over Row\" property of radar viewer to \"false\"", () => setProperty(page, "Show Mouse Over Row", el("radar viewer"), "false"));
      await session.step(26, "And user clicks on the \"cell 5 of AGE\" area of grid", () => clickArea(page, "cell 5 of AGE", el("grid")));
      await session.step(27, "Then row 5 should be current", () => currentRowIs(page, 5));
      await session.step(28, "And the \"view\" area of radar viewer should contain the color \"#00FF00\"", () => areaColor(page, "view", el("radar viewer"), "#00FF00"));
      await session.step(29, "When user takes a snapshot of radar viewer", () => takeSnapshot(page, el("radar viewer")));
      await session.step(30, "And user clicks on the \"cell 9 of AGE\" area of grid", () => clickArea(page, "cell 9 of AGE", el("grid")));
      await session.step(31, "Then row 9 should be current", () => currentRowIs(page, 9));
      await session.step(32, "And radar viewer should have repainted", () => repainted(page, el("radar viewer")));
      await session.step(33, "When user clicks on settings icon of radar viewer", () => clickOn(page, el("settings icon of radar viewer")));
      await session.step(34, "Given \"Selection\" category in context panel is expanded", () => isExpanded(page, el("\"Selection\" category in context panel")));
      await session.step(35, "When user unchecks \"Show Current Row\" property in context panel", () => uncheck(page, el("\"Show Current Row\" property in context panel")));
      await session.step(36, "Then \"Show Current Row\" property of radar viewer should be \"false\"", () => propertyShouldBe(page, "Show Current Row", el("radar viewer"), "false"));
      await session.step(37, "And the \"view\" area of radar viewer should not contain the color \"#00FF00\"", () => areaNotColor(page, "view", el("radar viewer"), "#00FF00"));
      await session.step(38, "When user checks \"Show Current Row\" property in context panel", () => check(page, el("\"Show Current Row\" property in context panel")));
      await session.step(39, "Then the \"view\" area of radar viewer should contain the color \"#00FF00\"", () => areaColor(page, "view", el("radar viewer"), "#00FF00"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pick Up / Apply copies a setting between two Radar viewers, and the two stay apart (public#1472)", async () => {
      await session.step(44, "Given user adds a radar viewer", () => addViewer(page, "radar"));
      await session.step(45, "Then the open tableview should have 2 radar viewers", () => viewerCount(page, 2, "radar"));
      await session.step(46, "When user sets \"normalization\" property of first radar viewer to \"Global\"", () => setProperty(page, "normalization", el("first radar viewer"), "Global"));
      await session.step(47, "Then \"normalization\" property of second radar viewer should be \"Column\"", () => propertyShouldBe(page, "normalization", el("second radar viewer"), "Column"));
      await session.step(48, "When user picks \"Pick Up / Apply > Pick Up\" from the viewer menu of first radar viewer", () => pickFromViewerMenu(page, "Pick Up / Apply > Pick Up", el("first radar viewer")));
      await session.step(49, "And user picks \"Pick Up / Apply > Apply\" from the viewer menu of second radar viewer", () => pickFromViewerMenu(page, "Pick Up / Apply > Apply", el("second radar viewer")));
      await session.step(50, "Then \"normalization\" property of second radar viewer should be \"Global\"", () => propertyShouldBe(page, "normalization", el("second radar viewer"), "Global"));
      await session.step(51, "When user takes a snapshot of first radar viewer", () => takeSnapshot(page, el("first radar viewer")));
      await session.step(52, "And user sets \"normalization\" property of first radar viewer to \"Column\"", () => setProperty(page, "normalization", el("first radar viewer"), "Column"));
      await session.step(53, "Then first radar viewer should have repainted", () => repainted(page, el("first radar viewer")));
      await session.step(54, "And \"normalization\" property of second radar viewer should be \"Global\"", () => propertyShouldBe(page, "normalization", el("second radar viewer"), "Global"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
