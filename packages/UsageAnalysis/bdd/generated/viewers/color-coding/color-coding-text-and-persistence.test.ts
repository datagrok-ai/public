/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/color-coding/color-coding-text-and-persistence.feature
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
import {closeCurrentView, layoutsDeleted, openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {applyNamedLayout, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingsDiffer, readingsEqual, saveNamedLayout} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A coloring applied to the text survives Extract and a saved layout, and the grid-wide All colours a boolean column", () => {
  const session = feature(test, "features/viewers/color-coding/color-coding-text-and-persistence.feature", import.meta.url);
  test("A coloring applied to the text survives Extract and a saved layout, and the grid-wide All colours a boolean column", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the layouts named \"bdd-cc-text-{run}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-cc-text-{run}")));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await run.scenario("Extract keeps the coloring on the text", async () => {
      await session.step(19, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(20, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(21, "And user selects \"text\" in Apply to input in \"Color-coding: AGE\" dialog", () => selectIn(page, "text", el("Apply to input in \"Color-coding: AGE\" dialog")));
      await session.step(22, "And user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(23, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(24, "And the text of \"AGE\" column should be color-coded", () => textColorCoded(page, "AGE"));
      await session.step(25, "When user picks \"Extract\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Extract", "header AGE", el("grid")));
      await session.step(26, "Then the open table views should be exactly \"demog-1000, result\"", () => openTableViewsExactly(page, "demog-1000, result"));
      await session.step(27, "When user switches to the \"result\" table view", () => switchTableView(page, "result"));
      await session.step(28, "Then the table should have the columns \"AGE\"", () => columnsExactly(page, "AGE"));
      await session.step(29, "And \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
      await session.step(30, "And the text of \"AGE\" column should be color-coded", () => textColorCoded(page, "AGE"));
      await session.step(31, "When user closes the current view", () => closeCurrentView(page));
      await session.step(32, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(33, "Then the open table views should be exactly \"demog-1000\"", () => openTableViewsExactly(page, "demog-1000"));
      await session.step(34, "When user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(35, "And user selects \"background\" in Apply to input in \"Color-coding: AGE\" dialog", () => selectIn(page, "background", el("Apply to input in \"Color-coding: AGE\" dialog")));
      await session.step(36, "And user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(37, "And user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(38, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved layout brings the text coloring back after it was switched off", async () => {
      await session.step(42, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(43, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(44, "And user selects \"text\" in Apply to input in \"Color-coding: AGE\" dialog", () => selectIn(page, "text", el("Apply to input in \"Color-coding: AGE\" dialog")));
      await session.step(45, "And user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(46, "Then the text of \"AGE\" column should be color-coded", () => textColorCoded(page, "AGE"));
      await session.step(47, "When user saves the layout of the current table view as \"bdd-cc-text-{run}\"", () => saveNamedLayout(page, session.text("bdd-cc-text-{run}")));
      await session.step(48, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(49, "And user selects \"background\" in Apply to input in \"Color-coding: AGE\" dialog", () => selectIn(page, "background", el("Apply to input in \"Color-coding: AGE\" dialog")));
      await session.step(50, "And user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(51, "And user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(52, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(53, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(54, "When user applies the layout \"bdd-cc-text-{run}\" to the current table view", () => applyNamedLayout(page, session.text("bdd-cc-text-{run}")));
      await session.step(55, "Then \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
      await session.step(56, "And the text of \"AGE\" column should be color-coded", () => textColorCoded(page, "AGE"));
      await session.step(57, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(58, "When user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(59, "And user selects \"background\" in Apply to input in \"Color-coding: AGE\" dialog", () => selectIn(page, "background", el("Apply to input in \"Color-coding: AGE\" dialog")));
      await session.step(60, "And user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(61, "And user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(62, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Grid Color Coding > All colours a boolean column too", async () => {
      await session.step(66, "Then \"CONTROL\" column should have no color coding", () => noColorCoding(page, "CONTROL"));
      await session.step(67, "And the \"color of cell 1 of CONTROL\" and \"color of cell 3 of CONTROL\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of CONTROL", "color of cell 3 of CONTROL", el("grid")));
      await session.step(68, "When user picks \"Grid Color Coding > All\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > All", "cell 1 of AGE", el("grid")));
      await session.step(69, "Then \"Color Coding\" property of grid should be \"All\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "All"));
      await session.step(70, "And the \"color of cell 1 of CONTROL\" and \"color of cell 3 of CONTROL\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of CONTROL", "color of cell 3 of CONTROL", el("grid")));
      await session.step(71, "When user picks \"Grid Color Coding > Auto\" from the context menu of the \"cell 1 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Grid Color Coding > Auto", "cell 1 of AGE", el("grid")));
      await session.step(72, "Then \"Color Coding\" property of grid should be \"Auto\"", () => propertyShouldBe(page, "Color Coding", el("grid"), "Auto"));
      await session.step(73, "And the \"color of cell 1 of CONTROL\" and \"color of cell 3 of CONTROL\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of CONTROL", "color of cell 3 of CONTROL", el("grid")));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
