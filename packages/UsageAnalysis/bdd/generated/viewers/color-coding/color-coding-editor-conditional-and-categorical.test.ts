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
import {clickOn, enterInto, hoverOver, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {categoricalColorIs, colorCodedAs, colorCodedCategorically, noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickColorSwatch, pickFromAreaContextMenu, readingNotAsRemembered, readingReads, readingsDiffer, readingsEqual, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The color coding editor of conditional and categorical columns: bins, rows, categories and Edit category colors", () => {
  const session = feature(test, "features/viewers/color-coding/color-coding-editor-conditional-and-categorical.feature", import.meta.url);
  test("The color coding editor of conditional and categorical columns: bins, rows, categories and Edit category colors", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await session.step(23, "And \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
    await run.scenario("The conditional editor opens on the default bins", async () => {
      await session.step(26, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(27, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(28, "When user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(29, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(30, "And Type input in \"Color-coding: AGE\" dialog should have value \"Conditional\"", () => shouldHaveValue(page, el("Type input in \"Color-coding: AGE\" dialog"), "Conditional"));
      await session.step(31, "And \"Bin...\" link in \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog"), "visible"));
      await session.step(32, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(33, "And the \"color of cell 1 of AGE\" and \"color of cell 4 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 4 of AGE", el("grid")));
      await session.step(34, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(35, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(36, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(37, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bin... re-bins the column into the number of bins asked for", async () => {
      await session.step(41, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(42, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(43, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(44, "When user clicks on \"Bin...\" link in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog")));
      await session.step(45, "Then \"Edit\" dialog should be visible", () => shouldBe(page, el("\"Edit\" dialog"), "visible"));
      await session.step(46, "When user enters \"2\" into \"Number of bins\" input in \"Edit\" dialog", () => enterInto(page, "2", el("\"Number of bins\" input in \"Edit\" dialog")));
      await session.step(47, "And user clicks on OK button in \"Edit\" dialog", () => clickOn(page, el("OK button in \"Edit\" dialog")));
      await session.step(48, "Then \"Edit\" dialog should be hidden", () => shouldBe(page, el("\"Edit\" dialog"), "hidden"));
      await session.step(49, "And \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(50, "And the \"color of cell 1 of AGE\" and \"color of cell 4 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 4 of AGE", el("grid")));
      await session.step(51, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(52, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(53, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(54, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(55, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Remove row takes a bin away", async () => {
      await session.step(59, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(60, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(61, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(62, "When user remembers the \"color of cell 3 of AGE\" reading of grid", () => rememberReading(page, "color of cell 3 of AGE", el("grid")));
      await session.step(63, "And user clicks on last Remove row button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("last Remove row button in \"Color-coding: AGE\" dialog")));
      await session.step(64, "Then the \"color of cell 3 of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(65, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(66, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(67, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(68, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A pattern typed into a row recolours the cells it now matches, and a bad pattern is marked", async () => {
      await session.step(72, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(73, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(74, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(75, "When user clicks on \"Bin...\" link in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog")));
      await session.step(76, "And user enters \"4\" into \"Number of bins\" input in \"Edit\" dialog", () => enterInto(page, "4", el("\"Number of bins\" input in \"Edit\" dialog")));
      await session.step(77, "And user clicks on OK button in \"Edit\" dialog", () => clickOn(page, el("OK button in \"Edit\" dialog")));
      await session.step(78, "Then \"Edit\" dialog should be hidden", () => shouldBe(page, el("\"Edit\" dialog"), "hidden"));
      await session.step(79, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(80, "When user enters \">20\" into \"Pattern 1\" text input in \"Color-coding: AGE\" dialog", () => enterInto(page, ">20", el("\"Pattern 1\" text input in \"Color-coding: AGE\" dialog")));
      await session.step(81, "Then the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(82, "When user enters \"abc\" into \"Pattern 1\" text input in \"Color-coding: AGE\" dialog", () => enterInto(page, "abc", el("\"Pattern 1\" text input in \"Color-coding: AGE\" dialog")));
      await session.step(83, "Then \"Pattern 1\" text input in \"Color-coding: AGE\" dialog should be invalid", () => shouldBe(page, el("\"Pattern 1\" text input in \"Color-coding: AGE\" dialog"), "invalid"));
      await session.step(84, "When user clicks on \"Pattern 2\" text input in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Pattern 2\" text input in \"Color-coding: AGE\" dialog")));
      await session.step(85, "And user hovers over \"Pattern 1\" text input in \"Color-coding: AGE\" dialog", () => hoverOver(page, el("\"Pattern 1\" text input in \"Color-coding: AGE\" dialog")));
      await session.step(86, "Then tooltip should contain text \"Invalid pattern\"", () => shouldContainText(page, el("tooltip"), "Invalid pattern"));
      await session.step(87, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(88, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(89, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(90, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A single value as a pattern colours the cells equal to it (public#1967)", async () => {
      await session.step(94, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(95, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(96, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(97, "When user clicks on \"Bin...\" link in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog")));
      await session.step(98, "And user enters \"4\" into \"Number of bins\" input in \"Edit\" dialog", () => enterInto(page, "4", el("\"Number of bins\" input in \"Edit\" dialog")));
      await session.step(99, "And user clicks on OK button in \"Edit\" dialog", () => clickOn(page, el("OK button in \"Edit\" dialog")));
      await session.step(100, "Then \"Edit\" dialog should be hidden", () => shouldBe(page, el("\"Edit\" dialog"), "hidden"));
      await session.step(101, "And the \"color of cell 1 of AGE\" and \"color of cell 2 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 2 of AGE", el("grid")));
      await session.step(102, "When user enters \"30\" into \"Pattern 1\" text input in \"Color-coding: AGE\" dialog", () => enterInto(page, "30", el("\"Pattern 1\" text input in \"Color-coding: AGE\" dialog")));
      await session.step(103, "Then the \"color of cell 1 of AGE\" and \"color of cell 2 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 2 of AGE", el("grid")));
      await session.step(104, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(105, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(106, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(107, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The categorical editor lists the categories, a colour picked for one is the colour the grid paints it, and Off clears the column", async () => {
      await session.step(111, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEX", el("grid")));
      await session.step(112, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(113, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(114, "When user picks \"Color Coding > Edit...\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header SEX", el("grid")));
      await session.step(115, "Then \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: SEX\" dialog"), "visible"));
      await session.step(116, "And \"F\" item in \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"F\" item in \"Color-coding: SEX\" dialog"), "visible"));
      await session.step(117, "And \"M\" item in \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"M\" item in \"Color-coding: SEX\" dialog"), "visible"));
      await session.step(118, "When user clicks on color bar of \"M\" item in \"Color-coding: SEX\" dialog", () => clickOn(page, el("color bar of \"M\" item in \"Color-coding: SEX\" dialog")));
      await session.step(119, "And user picks the color \"#9467BD\" in the color picker dialog", () => pickColorSwatch(page, "#9467BD"));
      await session.step(120, "Then the categorical color of \"M\" in \"SEX\" column should be \"#9467BD\"", () => categoricalColorIs(page, "M", "SEX", "#9467BD"));
      await session.step(121, "And the \"color of cell 4 of SEX\" reading of grid should be \"#9467bd\"", () => readingReads(page, "color of cell 4 of SEX", el("grid"), "#9467bd"));
      await session.step(122, "When user selects \"Off\" in Type input in \"Color-coding: SEX\" dialog", () => selectIn(page, "Off", el("Type input in \"Color-coding: SEX\" dialog")));
      await session.step(123, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(124, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(125, "When user clicks on CLOSE button in \"Color-coding: SEX\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: SEX\" dialog")));
      await session.step(126, "Then \"Color-coding: SEX\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: SEX\" dialog"), "hidden"));
      await session.step(127, "And \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(128, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Column Properties > Edit category colors... opens the categorical editor", async () => {
      await session.step(131, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEX", el("grid")));
      await session.step(132, "And user picks \"Column Properties...\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header SEX", el("grid")));
      await session.step(133, "Then \"SEX\" dialog should be visible", () => shouldBe(page, el("\"SEX\" dialog"), "visible"));
      await session.step(134, "When user clicks on \"Edit category colors...\" link in \"SEX\" dialog", () => clickOn(page, el("\"Edit category colors...\" link in \"SEX\" dialog")));
      await session.step(135, "Then \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: SEX\" dialog"), "visible"));
      await session.step(136, "And \"M\" item in \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"M\" item in \"Color-coding: SEX\" dialog"), "visible"));
      await session.step(137, "When user clicks on CLOSE button in \"Color-coding: SEX\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: SEX\" dialog")));
      await session.step(138, "And user clicks on OK button in \"SEX\" dialog", () => clickOn(page, el("OK button in \"SEX\" dialog")));
      await session.step(139, "Then \"SEX\" dialog should be hidden", () => shouldBe(page, el("\"SEX\" dialog"), "hidden"));
      await session.step(140, "When user picks \"Color Coding > Off\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header SEX", el("grid")));
      await session.step(141, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(142, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
