/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tooltips/line-chart-aggregated-tooltip.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tooltips]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {pickColumnInDialog} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {dialogCloses, openDataset, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hoverFirstArea, noBalloons, noErrors, pickFromContextMenu, pointerAway, propertyShouldBe, readingIs, readingReads, reportsNoError, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The line chart's aggregated tooltip with a split", () => {
  const session = feature(test, "features/viewers/tooltips/line-chart-aggregated-tooltip.feature", import.meta.url);
  test("The line chart's aggregated tooltip with a split", {tag: ["@journey", "@viewers", "@realizes:viewers.tooltips"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(19, "And user adds a line chart viewer", () => addViewer(page, "line chart"));
    await session.step(20, "And user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["X","Chemist 521"],["Y","CAST Idea ID"]]), [["X","Chemist 521"],["Y","CAST Idea ID"]]);
    await run.scenario("Edit Aggregated Tooltip adds two aggregations", async () => {
      await session.step(25, "Then the \"aggregated\" reading of line chart viewer should be \"true\"", () => readingReads(page, "aggregated", el("line chart viewer"), "true"));
      await session.step(26, "When user picks \"Tooltip > Edit...\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Tooltip > Edit...", el("line chart viewer")));
      await session.step(27, "Then \"Edit Aggregated Tooltip\" dialog should be visible", () => shouldBe(page, el("\"Edit Aggregated Tooltip\" dialog"), "visible"));
      await session.step(28, "When user clicks on first button in \"Edit Aggregated Tooltip\" dialog", () => clickOn(page, el("first button in \"Edit Aggregated Tooltip\" dialog")));
      await session.step(29, "And user picks \"Stereo Category\" in the column selector of the \"Edit Aggregated Tooltip\" dialog", () => pickColumnInDialog(page, "Stereo Category", "Edit Aggregated Tooltip"));
      await session.step(30, "And user selects \"concat unique\" in choice input in \"Edit Aggregated Tooltip\" dialog", () => selectIn(page, "concat unique", el("choice input in \"Edit Aggregated Tooltip\" dialog")));
      await session.step(31, "And user clicks on \"Add aggregation\" button in \"Edit Aggregated Tooltip\" dialog", () => clickOn(page, el("\"Add aggregation\" button in \"Edit Aggregated Tooltip\" dialog")));
      await session.step(32, "And user selects \"Average Mass\" in second column selector in \"Edit Aggregated Tooltip\" dialog", () => selectIn(page, "Average Mass", el("second column selector in \"Edit Aggregated Tooltip\" dialog")));
      await session.step(33, "And user selects \"min\" in second choice input in \"Edit Aggregated Tooltip\" dialog", () => selectIn(page, "min", el("second choice input in \"Edit Aggregated Tooltip\" dialog")));
      await session.step(34, "And user clicks on OK button in \"Edit Aggregated Tooltip\" dialog", () => clickOn(page, el("OK button in \"Edit Aggregated Tooltip\" dialog")));
      await session.step(35, "Then the \"Edit Aggregated Tooltip\" dialog should close", () => dialogCloses(page, "Edit Aggregated Tooltip"));
      await session.step(36, "And \"aggTooltipColumns\" property of line chart viewer should be \"concat unique(Stereo Category)\\nmin(Average Mass)\"", () => propertyShouldBe(page, "aggTooltipColumns", el("line chart viewer"), "concat unique(Stereo Category)\\nmin(Average Mass)"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Split by Stereo Category, a point's tooltip shows the two aggregations", async () => {
      await session.step(40, "When user sets \"Split\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "Split", el("line chart viewer"), "Stereo Category"));
      await session.step(41, "Then the \"split columns\" reading of line chart viewer should be 1", () => readingIs(page, "split columns", el("line chart viewer"), 1));
      await session.step(42, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(43, "And the \"markers drawn\" reading of line chart viewer should be 32", () => readingIs(page, "markers drawn", el("line chart viewer"), 32));
      await session.step(44, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(45, "When user hovers over the first \"point\" area of line chart viewer", () => hoverFirstArea(page, "point", el("line chart viewer")));
      await session.step(46, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(47, "And tooltip should contain text \"concat unique(Stereo Category)\"", () => shouldContainText(page, el("tooltip"), "concat unique(Stereo Category)"));
      await session.step(48, "And tooltip should contain text \"min(Average Mass)\"", () => shouldContainText(page, el("tooltip"), "min(Average Mass)"));
      await session.step(49, "When user moves the pointer away from line chart viewer", () => pointerAway(page, el("line chart viewer")));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
      await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
