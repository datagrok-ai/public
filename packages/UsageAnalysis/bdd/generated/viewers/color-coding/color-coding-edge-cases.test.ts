/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/color-coding/color-coding-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, core#1139, core#1171]
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
import {clickOn, enterInto, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {colorCodedAs, colorCodedCategorically, colorLinear, colorLinked, colorLinkedTo, columnsSelected, noColorCoding, removeColumn, renameColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickAreaHolding, closeContextMenu, noErrors, pickFromAreaContextMenu, readingsEqual, rightClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Color coding edge cases", () => {
  const session = feature(test, "features/viewers/color-coding/color-coding-edge-cases.feature", import.meta.url);
  test("A categorical coloring picked up cannot be applied to a numeric column", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await session.step(15, "And \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
    await session.step(19, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header RACE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header RACE", el("grid")));
    await session.step(20, "Then \"RACE\" column should be color-coded categorically", () => colorCodedCategorically(page, "RACE"));
    await session.step(21, "When user picks \"Color Coding > Pick Up Coloring\" from the context menu of the \"header RACE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Pick Up Coloring", "header RACE", el("grid")));
    await session.step(22, "And user right-clicks on the \"header SEX\" area of grid", () => rightClickArea(page, "header SEX", el("grid")));
    await session.step(23, "And user hovers over \"Color Coding\" menu item", () => hoverOver(page, el("\"Color Coding\" menu item")));
    await session.step(24, "Then \"Color-Coding---Apply-Coloring\" menu item should be enabled", () => shouldBe(page, el("\"Color-Coding---Apply-Coloring\" menu item"), "enabled"));
    await session.step(25, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(26, "And user right-clicks on the \"header AGE\" area of grid", () => rightClickArea(page, "header AGE", el("grid")));
    await session.step(27, "And user hovers over \"Color Coding\" menu item", () => hoverOver(page, el("\"Color Coding\" menu item")));
    await session.step(28, "Then \"Color-Coding---Apply-Coloring\" menu item should be disabled", () => shouldBe(page, el("\"Color-Coding---Apply-Coloring\" menu item"), "disabled"));
    await session.step(29, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(30, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
  });
  test("A linear coloring applied to a selection of columns colours only the numeric ones", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await session.step(15, "And \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
    await session.step(35, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
    await session.step(36, "Then \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
    await session.step(37, "When user picks \"Color Coding > Pick Up Coloring\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Pick Up Coloring", "header AGE", el("grid")));
    await session.step(38, "And user clicks on the \"header SEX\" area of grid holding Control", () => clickAreaHolding(page, "header SEX", el("grid"), "Control"));
    await session.step(39, "And user clicks on the \"header HEIGHT\" area of grid holding Control", () => clickAreaHolding(page, "header HEIGHT", el("grid"), "Control"));
    await session.step(40, "Then columns \"SEX, HEIGHT\" should be selected", () => columnsSelected(page, "SEX, HEIGHT"));
    await session.step(41, "When user picks \"Color Coding > Apply Coloring\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Apply Coloring", "header HEIGHT", el("grid")));
    await session.step(42, "Then \"HEIGHT\" column should be color-coded linearly", () => colorCodedAs(page, "HEIGHT", "linearly"));
    await session.step(43, "And \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
  });
  test("A linked coloring keeps following its source column after the source is renamed", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@realizes:core#1139"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await session.step(15, "And \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
    await session.step(48, "When user colors \"AGE\" column linearly from \"#0000FF\" to \"#FF0000\"", () => colorLinear(page, "AGE", "#0000FF", "#FF0000"));
    await session.step(49, "And user colors \"RACE\" column linked to \"AGE\" column", () => colorLinked(page, "RACE", "AGE"));
    await session.step(50, "Then the coloring of \"RACE\" column should be linked to \"AGE\" column", () => colorLinkedTo(page, "RACE", "AGE"));
    await session.step(51, "And the \"color of cell 1 of RACE\" and \"color of cell 1 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of RACE", "color of cell 1 of AGE", el("grid")));
    await session.step(52, "And the \"color of cell 2 of RACE\" and \"color of cell 2 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 2 of RACE", "color of cell 2 of AGE", el("grid")));
    await session.step(53, "When user renames \"AGE\" column to \"AGE2\"", () => renameColumn(page, "AGE", "AGE2"));
    await session.step(54, "Then \"AGE2\" column should be color-coded linearly", () => colorCodedAs(page, "AGE2", "linearly"));
    await session.step(55, "And the \"color of cell 2 of RACE\" and \"color of cell 2 of AGE2\" readings of grid should be the same", () => readingsEqual(page, "color of cell 2 of RACE", "color of cell 2 of AGE2", el("grid")));
  });
  test("Bin... refuses a number of bins below one and keeps the coloring", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:66", "@realizes:core#1171"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await session.step(15, "And \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
    await knownFailure(async () => {
      await session.step(59, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(60, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(61, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(62, "When user clicks on \"Bin...\" link in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog")));
      await session.step(63, "Then \"Edit\" dialog should be visible", () => shouldBe(page, el("\"Edit\" dialog"), "visible"));
      await session.step(64, "When user enters \"0\" into \"Number of bins\" input in \"Edit\" dialog", () => enterInto(page, "0", el("\"Number of bins\" input in \"Edit\" dialog")));
      await session.step(65, "And user clicks on OK button in \"Edit\" dialog", () => clickOn(page, el("OK button in \"Edit\" dialog")));
      await session.step(66, "Then \"Edit\" dialog should be visible", () => shouldBe(page, el("\"Edit\" dialog"), "visible"));
      await session.step(67, "When user clicks on CANCEL button in \"Edit\" dialog", () => clickOn(page, el("CANCEL button in \"Edit\" dialog")));
      await session.step(68, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    }, 66);
  });
  test("A linked coloring whose source column is removed paints nothing", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@realizes:core#1139"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await session.step(15, "And \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
    await session.step(73, "When user colors \"AGE\" column linearly from \"#0000FF\" to \"#FF0000\"", () => colorLinear(page, "AGE", "#0000FF", "#FF0000"));
    await session.step(74, "And user colors \"RACE\" column linked to \"AGE\" column", () => colorLinked(page, "RACE", "AGE"));
    await session.step(75, "Then the \"color of cell 1 of RACE\" and \"color of cell 1 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of RACE", "color of cell 1 of AGE", el("grid")));
    await session.step(76, "When user removes \"AGE\" column", () => removeColumn(page, "AGE"));
    await session.step(77, "Then the \"color of cell 1 of RACE\" and \"color of cell 1 of SEX\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of RACE", "color of cell 1 of SEX", el("grid")));
    await session.step(78, "And the \"color of cell 2 of RACE\" and \"color of cell 2 of SEX\" readings of grid should be the same", () => readingsEqual(page, "color of cell 2 of RACE", "color of cell 2 of SEX", el("grid")));
    await session.step(79, "And no errors should have been logged", () => noErrors(page));
  });
});
