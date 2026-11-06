/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-category-order-molecules.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
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
import {pickColumnInDialog} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset, packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noBalloons, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLiesBeside, areasLeftToRightOnceComputed} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot molecule categories ordered by a computed property", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-category-order-molecules.feature", import.meta.url);
  test("Scatter plot molecule categories ordered by a computed property", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","Core"],["Y","Series"]]), [["X","Core"],["Y","Series"]]);
    await session.step(22, "Then the \"marker of row 2\" area of scatter plot viewer should lie to the left of the \"marker of row 1\" area", () => areaLiesBeside(page, "marker of row 2", el("scatter plot viewer"), "left", "marker of row 1"));
    await session.step(23, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 8\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 8"));
    await run.scenario("Molecular weight orders the cores, in both directions", async () => {
      await session.step(26, "When user picks \"Order...\" from the context menu of the \"left edge of x axis\" area of scatter plot viewer", () => pickFromAreaContextMenu(page, "Order...", "left edge of x axis", el("scatter plot viewer")));
      await session.step(27, "Then \"Order categories: Core\" dialog should be visible", () => shouldBe(page, el("\"Order categories: Core\" dialog"), "visible"));
      await session.step(28, "When user selects \"Aggregated value\" in \"Order by\" input in \"Order categories: Core\" dialog", () => selectIn(page, "Aggregated value", el("\"Order by\" input in \"Order categories: Core\" dialog")));
      await session.step(29, "Then \"Column\" input in \"Order categories: Core\" dialog should contain text \"Structure\"", () => shouldContainText(page, el("\"Column\" input in \"Order categories: Core\" dialog"), "Structure"));
      await session.step(30, "When user picks \"Core\" in the column selector of the \"Order categories: Core\" dialog", () => pickColumnInDialog(page, "Core", "Order categories: Core"));
      await session.step(31, "Then \"Property\" input in \"Order categories: Core\" dialog should be visible", () => shouldBe(page, el("\"Property\" input in \"Order categories: Core\" dialog"), "visible"));
      await session.step(32, "And \"Property\" input in \"Order categories: Core\" dialog should have value \"MW\"", () => shouldHaveValue(page, el("\"Property\" input in \"Order categories: Core\" dialog"), "MW"));
      await session.step(33, "And \"Aggregation\" input in \"Order categories: Core\" dialog should be hidden", () => shouldBe(page, el("\"Aggregation\" input in \"Order categories: Core\" dialog"), "hidden"));
      await session.step(34, "When user selects \"Ascending\" in \"Direction\" input in \"Order categories: Core\" dialog", () => selectIn(page, "Ascending", el("\"Direction\" input in \"Order categories: Core\" dialog")));
      await session.step(35, "Then once the computed order arrives, the areas \"marker of row 8, marker of row 1, marker of row 2\" of scatter plot viewer should lie from left to right", () => areasLeftToRightOnceComputed(page, "marker of row 8, marker of row 1, marker of row 2", el("scatter plot viewer")));
      await session.step(36, "When user selects \"Descending\" in \"Direction\" input in \"Order categories: Core\" dialog", () => selectIn(page, "Descending", el("\"Direction\" input in \"Order categories: Core\" dialog")));
      await session.step(37, "Then the \"marker of row 2\" area of scatter plot viewer should lie to the left of the \"marker of row 1\" area", () => areaLiesBeside(page, "marker of row 2", el("scatter plot viewer"), "left", "marker of row 1"));
      await session.step(38, "And the \"marker of row 1\" area of scatter plot viewer should lie to the left of the \"marker of row 8\" area", () => areaLiesBeside(page, "marker of row 1", el("scatter plot viewer"), "left", "marker of row 8"));
      await session.step(39, "When user selects \"Default\" in \"Order by\" input in \"Order categories: Core\" dialog", () => selectIn(page, "Default", el("\"Order by\" input in \"Order categories: Core\" dialog")));
      await session.step(40, "And user clicks on \"CLOSE\" button in \"Order categories: Core\" dialog", () => clickOn(page, el("\"CLOSE\" button in \"Order categories: Core\" dialog")));
      await session.step(41, "Then the \"marker of row 2\" area of scatter plot viewer should lie to the left of the \"marker of row 1\" area", () => areaLiesBeside(page, "marker of row 2", el("scatter plot viewer"), "left", "marker of row 1"));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
