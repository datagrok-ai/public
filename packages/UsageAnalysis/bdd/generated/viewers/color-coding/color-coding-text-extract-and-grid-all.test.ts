/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/color-coding/color-coding-text-extract-and-grid-all.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
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
import {clickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorCodedAs, noColorCoding, textColorCoded} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors, pickFromAreaContextMenu, propertyShouldBe, readingsDiffer, readingsEqual} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A coloring applied to the text survives Extract, and the grid-wide All colours a boolean column", () => {
  const session = feature(test, "features/viewers/color-coding/color-coding-text-extract-and-grid-all.feature", import.meta.url);
  test("A coloring applied to the text survives Extract, and the grid-wide All colours a boolean column", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await run.scenario("Extract keeps the coloring on the text", async () => {
      await session.step(17, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(18, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(19, "And user selects \"text\" in Apply to input in \"Color-coding: AGE\" dialog", () => selectIn(page, "text", el("Apply to input in \"Color-coding: AGE\" dialog")));
      await session.step(20, "And user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(21, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(22, "And the text of \"AGE\" column should be color-coded", () => textColorCoded(page, "AGE"));
      await session.step(23, "When user picks \"Extract\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Extract", "header AGE", el("grid")));
      await session.step(24, "Then the open table views should be exactly \"demog-1000, result\"", () => openTableViewsExactly(page, "demog-1000, result"));
      await session.step(25, "When user switches to the \"result\" table view", () => switchTableView(page, "result"));
      await session.step(26, "Then the table should have the columns \"AGE\"", () => columnsExactly(page, "AGE"));
      await session.step(27, "And \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
      await session.step(28, "And the text of \"AGE\" column should be color-coded", () => textColorCoded(page, "AGE"));
      await session.step(29, "When user closes the current view", () => closeCurrentView(page));
      await session.step(30, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(31, "Then the open table views should be exactly \"demog-1000\"", () => openTableViewsExactly(page, "demog-1000"));
      await session.step(32, "When user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(33, "And user selects \"background\" in Apply to input in \"Color-coding: AGE\" dialog", () => selectIn(page, "background", el("Apply to input in \"Color-coding: AGE\" dialog")));
      await session.step(34, "And user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(35, "And user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(36, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Grid Color Coding > All colours a boolean column too", async () => {
      await session.step(40, "Then \"CONTROL\" column should have no color coding", () => noColorCoding(page, "CONTROL"));
      await session.step(41, "And the \"color of cell 1 of CONTROL\" and \"color of cell 3 of CONTROL\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of CONTROL", "color of cell 3 of CONTROL", el("grid")));
      await session.step(42, "When user picks \"Grid Color Coding > All\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > All", "cell 1 of AGE", el("grid")));
      await session.step(43, "Then \"Color Coding\" property of grid should be \"All\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "All"));
      await session.step(44, "And the \"color of cell 1 of CONTROL\" and \"color of cell 3 of CONTROL\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of CONTROL", "color of cell 3 of CONTROL", el("grid")));
      await session.step(45, "When user picks \"Grid Color Coding > Auto\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > Auto", "cell 1 of AGE", el("grid")));
      await session.step(46, "Then \"Color Coding\" property of grid should be \"Auto\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "Auto"));
      await session.step(47, "And the \"color of cell 1 of CONTROL\" and \"color of cell 3 of CONTROL\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of CONTROL", "color of cell 3 of CONTROL", el("grid")));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
