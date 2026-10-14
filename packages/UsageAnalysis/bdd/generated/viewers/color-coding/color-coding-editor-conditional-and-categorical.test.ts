/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/color-coding/color-coding-editor-conditional-and-categorical.feature
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
import {clickOn, enterInto, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {colorCodedAs, colorCodedCategorically, noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromAreaContextMenu, readingNotAsRemembered, readingsDiffer, readingsEqual, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The color coding editor of conditional and categorical columns: bins, rows, categories and Edit category colors", () => {
  const session = feature(test, "features/viewers/color-coding/color-coding-editor-conditional-and-categorical.feature", import.meta.url);
  test("The color coding editor of conditional and categorical columns: bins, rows, categories and Edit category colors", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await session.step(21, "And \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
    await run.scenario("The conditional editor opens on the default bins", async () => {
      await session.step(24, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(25, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(26, "When user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(27, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(28, "And Type input in \"Color-coding: AGE\" dialog should have value \"Conditional\"", () => shouldHaveValue(page, el("Type input in \"Color-coding: AGE\" dialog"), "Conditional"));
      await session.step(29, "And \"Bin...\" link in \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog"), "visible"));
      await session.step(30, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(31, "And the \"color of cell 1 of AGE\" and \"color of cell 4 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 4 of AGE", el("grid")));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bin... re-bins the column into the number of bins asked for", async () => {
      await session.step(35, "When user clicks on \"Bin...\" link in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog")));
      await session.step(36, "Then \"Edit\" dialog should be visible", () => shouldBe(page, el("\"Edit\" dialog"), "visible"));
      await session.step(37, "When user enters \"2\" into \"Number of bins\" input in \"Edit\" dialog", () => enterInto(page, "2", el("\"Number of bins\" input in \"Edit\" dialog")));
      await session.step(38, "And user clicks on OK button in \"Edit\" dialog", () => clickOn(page, el("OK button in \"Edit\" dialog")));
      await session.step(39, "Then \"Edit\" dialog should be hidden", () => shouldBe(page, el("\"Edit\" dialog"), "hidden"));
      await session.step(40, "And \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(41, "And the \"color of cell 1 of AGE\" and \"color of cell 4 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 4 of AGE", el("grid")));
      await session.step(42, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Remove row takes a bin away", async () => {
      await session.step(46, "When user remembers the \"color of cell 3 of AGE\" reading of grid", () => rememberReading(page, "color of cell 3 of AGE", el("grid")));
      await session.step(47, "And user clicks on last Remove row button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("last Remove row button in \"Color-coding: AGE\" dialog")));
      await session.step(48, "Then the \"color of cell 3 of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(49, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(50, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(51, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(52, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The categorical editor lists the categories, and Off in it clears the column", async () => {
      await session.step(56, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEX", el("grid")));
      await session.step(57, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(58, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(59, "When user picks \"Color Coding > Edit...\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header SEX", el("grid")));
      await session.step(60, "Then \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: SEX\" dialog"), "visible"));
      await session.step(61, "And \"F\" item in \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"F\" item in \"Color-coding: SEX\" dialog"), "visible"));
      await session.step(62, "And \"M\" item in \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"M\" item in \"Color-coding: SEX\" dialog"), "visible"));
      await session.step(63, "When user selects \"Off\" in Type input in \"Color-coding: SEX\" dialog", () => selectIn(page, "Off", el("Type input in \"Color-coding: SEX\" dialog")));
      await session.step(64, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(65, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(66, "When user clicks on CLOSE button in \"Color-coding: SEX\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: SEX\" dialog")));
      await session.step(67, "Then \"Color-coding: SEX\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: SEX\" dialog"), "hidden"));
      await session.step(68, "And \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Column Properties > Edit category colors... opens the categorical editor", async () => {
      await session.step(72, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEX", el("grid")));
      await session.step(73, "And user picks \"Column Properties...\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header SEX", el("grid")));
      await session.step(74, "Then \"SEX\" dialog should be visible", () => shouldBe(page, el("\"SEX\" dialog"), "visible"));
      await session.step(75, "When user clicks on \"Edit category colors...\" link in \"SEX\" dialog", () => clickOn(page, el("\"Edit category colors...\" link in \"SEX\" dialog")));
      await session.step(76, "Then \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: SEX\" dialog"), "visible"));
      await session.step(77, "And \"M\" item in \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"M\" item in \"Color-coding: SEX\" dialog"), "visible"));
      await session.step(78, "When user clicks on CLOSE button in \"Color-coding: SEX\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: SEX\" dialog")));
      await session.step(79, "And user clicks on OK button in \"SEX\" dialog", () => clickOn(page, el("OK button in \"SEX\" dialog")));
      await session.step(80, "Then \"SEX\" dialog should be hidden", () => shouldBe(page, el("\"SEX\" dialog"), "hidden"));
      await session.step(81, "When user picks \"Color Coding > Off\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header SEX", el("grid")));
      await session.step(82, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
