/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/density-plot/density-plot-bin-selector-escape-tooltip.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.density-plot]
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
import {pressKey, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, noneSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hasArea, hoverArea, noErrors, pointerAway, readingIs, readingReads, readingsEqual, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Density plot bin selector, Escape and the bin tooltip", () => {
  const session = feature(test, "features/viewers/density-plot/density-plot-bin-selector-escape-tooltip.feature", import.meta.url);
  test("Density plot bin selector, Escape and the bin tooltip", {tag: ["@journey", "@viewers", "@realizes:viewers.density-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a density plot viewer with:", () => addViewerWith(page, "density plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"]]);
    await session.step(16, "Then the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
    await session.step(17, "And density plot viewer should have a \"densest bin\" area", () => hasArea(page, el("density plot viewer"), "densest bin"));
    await run.scenario("Show Bin Selector takes the bin selector out of the layout and puts it back, as the layout resolved it", async () => {
      await session.step(20, "Then the \"bin selector shown\" reading of density plot viewer should be \"true\"", () => readingReads(page, "bin selector shown", el("density plot viewer"), "true"));
      await session.step(21, "When user sets \"showBinSelector\" property of density plot viewer to \"false\"", () => setProperty(page, "showBinSelector", el("density plot viewer"), "false"));
      await session.step(22, "Then the \"bin selector shown\" reading of density plot viewer should be \"false\"", () => readingReads(page, "bin selector shown", el("density plot viewer"), "false"));
      await session.step(23, "And the \"x selector shown\" reading of density plot viewer should be \"true\"", () => readingReads(page, "x selector shown", el("density plot viewer"), "true"));
      await session.step(24, "When user sets \"showBinSelector\" property of density plot viewer to \"true\"", () => setProperty(page, "showBinSelector", el("density plot viewer"), "true"));
      await session.step(25, "Then the \"bin selector shown\" reading of density plot viewer should be \"true\"", () => readingReads(page, "bin selector shown", el("density plot viewer"), "true"));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Escape clears the selection a bin click made", async () => {
      await session.step(29, "Given user clears the row selection", () => clearSelection(page));
      await session.step(30, "When user clicks on the \"densest bin\" area of density plot viewer", () => clickArea(page, "densest bin", el("density plot viewer")));
      await session.step(31, "Then the \"rows selected\" and \"rows in densest bin\" readings of density plot viewer should be the same", () => readingsEqual(page, "rows selected", "rows in densest bin", el("density plot viewer")));
      await session.step(32, "And some rows should be selected", () => someSelected(page));
      await session.step(33, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(34, "Then no rows should be selected", () => noneSelected(page));
      await session.step(35, "And the \"rows selected\" reading of density plot viewer should be 0", () => readingIs(page, "rows selected", el("density plot viewer"), 0));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A bin's tooltip gives the bin's range of both binned columns and its row count", async () => {
      await session.step(39, "Then the \"rows in densest bin\" reading of density plot viewer should be 5", () => readingIs(page, "rows in densest bin", el("density plot viewer"), 5));
      await session.step(40, "When user hovers over the \"densest bin\" area of density plot viewer", () => hoverArea(page, "densest bin", el("density plot viewer")));
      await session.step(41, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(42, "And tooltip should contain text \"AGE: 47..47\"", () => shouldContainText(page, el("tooltip"), "AGE: 47..47"));
      await session.step(43, "And tooltip should contain text \"WEIGHT: 62.20..64.40\"", () => shouldContainText(page, el("tooltip"), "WEIGHT: 62.20..64.40"));
      await session.step(44, "And tooltip should contain text \"5 rows\"", () => shouldContainText(page, el("tooltip"), "5 rows"));
      await session.step(45, "When user moves the pointer away from density plot viewer", () => pointerAway(page, el("density plot viewer")));
      await session.step(46, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
