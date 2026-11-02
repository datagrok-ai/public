/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/color-coding/color-coding-editor-conditional-and-categorical.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, core#1171]
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
import {closeCurrentView, openDataset, openTableOf, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickColorSwatch, pickFromAreaContextMenu, readingNotAsRemembered, readingReads, readingsDiffer, readingsEqual, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The color coding editor of conditional and categorical columns: bins, rows, categories and Edit category colors", () => {
  const session = feature(test, "features/viewers/color-coding/color-coding-editor-conditional-and-categorical.feature", import.meta.url);
  test("The color coding editor of conditional and categorical columns: bins, rows, categories and Edit category colors", {tag: ["@journey", "@viewers", "@realizes:viewers.grid", "@edge-cases", "@realizes:core#1171"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(26, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await session.step(27, "And \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
    await run.scenario("The conditional editor opens on the default bins", async () => {
      await session.step(30, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(31, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(32, "When user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(33, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(34, "And Type input in \"Color-coding: AGE\" dialog should have value \"Conditional\"", () => shouldHaveValue(page, el("Type input in \"Color-coding: AGE\" dialog"), "Conditional"));
      await session.step(35, "And \"Bin...\" link in \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog"), "visible"));
      await session.step(36, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(37, "And the \"color of cell 1 of AGE\" and \"color of cell 4 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 4 of AGE", el("grid")));
      await session.step(38, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(39, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(40, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(41, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bin... re-bins the column into the number of bins asked for", async () => {
      await session.step(45, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(46, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(47, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(48, "When user clicks on \"Bin...\" link in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog")));
      await session.step(49, "Then \"Edit\" dialog should be visible", () => shouldBe(page, el("\"Edit\" dialog"), "visible"));
      await session.step(50, "When user enters \"2\" into \"Number of bins\" input in \"Edit\" dialog", () => enterInto(page, "2", el("\"Number of bins\" input in \"Edit\" dialog")));
      await session.step(51, "And user clicks on OK button in \"Edit\" dialog", () => clickOn(page, el("OK button in \"Edit\" dialog")));
      await session.step(52, "Then \"Edit\" dialog should be hidden", () => shouldBe(page, el("\"Edit\" dialog"), "hidden"));
      await session.step(53, "And \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(54, "And the \"color of cell 1 of AGE\" and \"color of cell 4 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 4 of AGE", el("grid")));
      await session.step(55, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(56, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(57, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(58, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(59, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Remove row takes a bin away", async () => {
      await session.step(63, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(64, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(65, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(66, "When user remembers the \"color of cell 3 of AGE\" reading of grid", () => rememberReading(page, "color of cell 3 of AGE", el("grid")));
      await session.step(67, "And user clicks on last Remove row button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("last Remove row button in \"Color-coding: AGE\" dialog")));
      await session.step(68, "Then the \"color of cell 3 of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(69, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(70, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(71, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(72, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A pattern typed into a row recolours the cells it now matches, and a bad pattern is marked", async () => {
      await session.step(76, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(77, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(78, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(79, "When user clicks on \"Bin...\" link in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog")));
      await session.step(80, "And user enters \"4\" into \"Number of bins\" input in \"Edit\" dialog", () => enterInto(page, "4", el("\"Number of bins\" input in \"Edit\" dialog")));
      await session.step(81, "And user clicks on OK button in \"Edit\" dialog", () => clickOn(page, el("OK button in \"Edit\" dialog")));
      await session.step(82, "Then \"Edit\" dialog should be hidden", () => shouldBe(page, el("\"Edit\" dialog"), "hidden"));
      await session.step(83, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(84, "When user enters \">20\" into \"Pattern 1\" text input in \"Color-coding: AGE\" dialog", () => enterInto(page, ">20", el("\"Pattern 1\" text input in \"Color-coding: AGE\" dialog")));
      await session.step(85, "Then the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(86, "When user enters \"abc\" into \"Pattern 1\" text input in \"Color-coding: AGE\" dialog", () => enterInto(page, "abc", el("\"Pattern 1\" text input in \"Color-coding: AGE\" dialog")));
      await session.step(87, "Then \"Pattern 1\" text input in \"Color-coding: AGE\" dialog should be invalid", () => shouldBe(page, el("\"Pattern 1\" text input in \"Color-coding: AGE\" dialog"), "invalid"));
      await session.step(88, "When user clicks on \"Pattern 2\" text input in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Pattern 2\" text input in \"Color-coding: AGE\" dialog")));
      await session.step(89, "And user hovers over \"Pattern 1\" text input in \"Color-coding: AGE\" dialog", () => hoverOver(page, el("\"Pattern 1\" text input in \"Color-coding: AGE\" dialog")));
      await session.step(90, "Then tooltip should contain text \"Invalid pattern\"", () => shouldContainText(page, el("tooltip"), "Invalid pattern"));
      await session.step(91, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(92, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(93, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(94, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A single value as a pattern colours the cells equal to it (public#1967)", async () => {
      await session.step(98, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(99, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(100, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(101, "When user clicks on \"Bin...\" link in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog")));
      await session.step(102, "And user enters \"4\" into \"Number of bins\" input in \"Edit\" dialog", () => enterInto(page, "4", el("\"Number of bins\" input in \"Edit\" dialog")));
      await session.step(103, "And user clicks on OK button in \"Edit\" dialog", () => clickOn(page, el("OK button in \"Edit\" dialog")));
      await session.step(104, "Then \"Edit\" dialog should be hidden", () => shouldBe(page, el("\"Edit\" dialog"), "hidden"));
      await session.step(105, "And the \"color of cell 1 of AGE\" and \"color of cell 2 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 2 of AGE", el("grid")));
      await session.step(106, "When user enters \"30\" into \"Pattern 1\" text input in \"Color-coding: AGE\" dialog", () => enterInto(page, "30", el("\"Pattern 1\" text input in \"Color-coding: AGE\" dialog")));
      await session.step(107, "Then the \"color of cell 1 of AGE\" and \"color of cell 2 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 2 of AGE", el("grid")));
      await session.step(108, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(109, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(110, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(111, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The categorical editor lists the categories, a colour picked for one is the colour the grid paints it, and Off clears the column", async () => {
      await session.step(115, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEX", el("grid")));
      await session.step(116, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(117, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(118, "When user picks \"Color Coding > Edit...\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header SEX", el("grid")));
      await session.step(119, "Then \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: SEX\" dialog"), "visible"));
      await session.step(120, "And \"F\" item in \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"F\" item in \"Color-coding: SEX\" dialog"), "visible"));
      await session.step(121, "And \"M\" item in \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"M\" item in \"Color-coding: SEX\" dialog"), "visible"));
      await session.step(122, "When user clicks on color bar of \"M\" item in \"Color-coding: SEX\" dialog", () => clickOn(page, el("color bar of \"M\" item in \"Color-coding: SEX\" dialog")));
      await session.step(123, "And user picks the color \"#9467BD\" in the color picker dialog", () => pickColorSwatch(page, "#9467BD"));
      await session.step(124, "Then the categorical color of \"M\" in \"SEX\" column should be \"#9467BD\"", () => categoricalColorIs(page, "M", "SEX", "#9467BD"));
      await session.step(125, "And the \"color of cell 4 of SEX\" reading of grid should be \"#9467bd\"", () => readingReads(page, "color of cell 4 of SEX", el("grid"), "#9467bd"));
      await session.step(126, "When user selects \"Off\" in Type input in \"Color-coding: SEX\" dialog", () => selectIn(page, "Off", el("Type input in \"Color-coding: SEX\" dialog")));
      await session.step(127, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(128, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(129, "When user clicks on CLOSE button in \"Color-coding: SEX\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: SEX\" dialog")));
      await session.step(130, "Then \"Color-coding: SEX\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: SEX\" dialog"), "hidden"));
      await session.step(131, "And \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Column Properties > Edit category colors... opens the categorical editor", async () => {
      await session.step(135, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEX", el("grid")));
      await session.step(136, "And user picks \"Column Properties...\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header SEX", el("grid")));
      await session.step(137, "Then \"SEX\" dialog should be visible", () => shouldBe(page, el("\"SEX\" dialog"), "visible"));
      await session.step(138, "When user clicks on \"Edit category colors...\" link in \"SEX\" dialog", () => clickOn(page, el("\"Edit category colors...\" link in \"SEX\" dialog")));
      await session.step(139, "Then \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: SEX\" dialog"), "visible"));
      await session.step(140, "And \"M\" item in \"Color-coding: SEX\" dialog should be visible", () => shouldBe(page, el("\"M\" item in \"Color-coding: SEX\" dialog"), "visible"));
      await session.step(141, "When user clicks on CLOSE button in \"Color-coding: SEX\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: SEX\" dialog")));
      await session.step(142, "And user clicks on OK button in \"SEX\" dialog", () => clickOn(page, el("OK button in \"SEX\" dialog")));
      await session.step(143, "Then \"SEX\" dialog should be hidden", () => shouldBe(page, el("\"SEX\" dialog"), "hidden"));
      await session.step(144, "When user picks \"Color Coding > Off\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header SEX", el("grid")));
      await session.step(145, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(146, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bin... refuses fewer than one bin, and CANCEL leaves the coloring as it was (core#1171)", async () => {
      await session.step(150, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(151, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(152, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(153, "When user clicks on \"Bin...\" link in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog")));
      await session.step(154, "Then \"Edit\" dialog should be visible", () => shouldBe(page, el("\"Edit\" dialog"), "visible"));
      await session.step(155, "When user enters \"0\" into \"Number of bins\" input in \"Edit\" dialog", () => enterInto(page, "0", el("\"Number of bins\" input in \"Edit\" dialog")));
      await session.step(156, "Then OK button in \"Edit\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Edit\" dialog"), "disabled"));
      await session.step(157, "When user clicks on CANCEL button in \"Edit\" dialog", () => clickOn(page, el("CANCEL button in \"Edit\" dialog")));
      await session.step(158, "Then \"Edit\" dialog should be hidden", () => shouldBe(page, el("\"Edit\" dialog"), "hidden"));
      await session.step(159, "And \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(160, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(161, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(162, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(163, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(164, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(165, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A value on the edge two Bin... ranges share goes to the lower bin", async () => {
      await session.step(169, "Given user opens a table \"bin-edge\" with:", () => openTableOf(page, "bin-edge", [["Name","Level"],["a","0"],["b","5"],["c","10"]]), [["Name","Level"],["a","0"],["b","5"],["c","10"]]);
      await session.step(174, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header Level\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header Level", el("grid")));
      await session.step(175, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header Level\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header Level", el("grid")));
      await session.step(176, "Then \"Color-coding: Level\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: Level\" dialog"), "visible"));
      await session.step(177, "When user clicks on \"Bin...\" link in \"Color-coding: Level\" dialog", () => clickOn(page, el("\"Bin...\" link in \"Color-coding: Level\" dialog")));
      await session.step(178, "Then \"Edit\" dialog should be visible", () => shouldBe(page, el("\"Edit\" dialog"), "visible"));
      await session.step(179, "When user enters \"2\" into \"Number of bins\" input in \"Edit\" dialog", () => enterInto(page, "2", el("\"Number of bins\" input in \"Edit\" dialog")));
      await session.step(180, "And user clicks on OK button in \"Edit\" dialog", () => clickOn(page, el("OK button in \"Edit\" dialog")));
      await session.step(181, "Then \"Level\" column should be color-coded conditionally", () => colorCodedAs(page, "Level", "conditionally"));
      await session.step(182, "And the \"color of cell 1 of Level\" and \"color of cell 3 of Level\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of Level", "color of cell 3 of Level", el("grid")));
      await session.step(183, "And the \"color of cell 1 of Level\" and \"color of cell 2 of Level\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of Level", "color of cell 2 of Level", el("grid")));
      await session.step(184, "When user clicks on CLOSE button in \"Color-coding: Level\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: Level\" dialog")));
      await session.step(185, "Then \"Color-coding: Level\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: Level\" dialog"), "hidden"));
      await session.step(186, "When user closes the current view", () => closeCurrentView(page));
      await session.step(187, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(188, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
