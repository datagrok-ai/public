/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/color-coding/color-coding-editor-linear.feature
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
import {check, clearField, clickOn, enterInto, hoverOver, selectIn, shouldBe, shouldHaveValue, shouldNotBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorCodedAs, colorSchemeIs, noColorCoding, textColorCoded} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, noErrors, pickFromAreaContextMenu, readingAsRemembered, readingNotAsRemembered, readingsDiffer, readingsEqual, rememberReading, rightClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The color coding editor of a numeric column: type, Apply to, Min and Max, invert and Edit scheme", () => {
  const session = feature(test, "features/viewers/color-coding/color-coding-editor-linear.feature", import.meta.url);
  test("The color coding editor of a numeric column: type, Apply to, Min and Max, invert and Edit scheme", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await run.scenario("Each Color Coding item of the header menu replaces the column's coding, and Off removes it", async () => {
      await session.step(24, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(25, "Then \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
      await session.step(26, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(27, "When user right-clicks on the \"header AGE\" area of grid", () => rightClickArea(page, "header AGE", el("grid")));
      await session.step(28, "And user hovers over \"Color Coding\" menu item", () => hoverOver(page, el("\"Color Coding\" menu item")));
      await session.step(29, "Then \"Color-Coding---Linear\" menu item should be selected", () => shouldBe(page, el("\"Color-Coding---Linear\" menu item"), "selected"));
      await session.step(30, "And \"Color-Coding---Conditional\" menu item should not be selected", () => shouldNotBe(page, el("\"Color-Coding---Conditional\" menu item"), "selected"));
      await session.step(31, "And \"Color-Coding---Off\" menu item should not be selected", () => shouldNotBe(page, el("\"Color-Coding---Off\" menu item"), "selected"));
      await session.step(32, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(33, "And user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(34, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(35, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(36, "When user right-clicks on the \"header AGE\" area of grid", () => rightClickArea(page, "header AGE", el("grid")));
      await session.step(37, "And user hovers over \"Color Coding\" menu item", () => hoverOver(page, el("\"Color Coding\" menu item")));
      await session.step(38, "Then \"Color-Coding---Conditional\" menu item should be selected", () => shouldBe(page, el("\"Color-Coding---Conditional\" menu item"), "selected"));
      await session.step(39, "And \"Color-Coding---Linear\" menu item should not be selected", () => shouldNotBe(page, el("\"Color-Coding---Linear\" menu item"), "selected"));
      await session.step(40, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(41, "And user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(42, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(43, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Type choice of the editor switches the column's coding", async () => {
      await session.step(47, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(48, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(49, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(50, "And Type input in \"Color-coding: AGE\" dialog should have value \"Linear\"", () => shouldHaveValue(page, el("Type input in \"Color-coding: AGE\" dialog"), "Linear"));
      await session.step(51, "When user selects \"Conditional\" in Type input in \"Color-coding: AGE\" dialog", () => selectIn(page, "Conditional", el("Type input in \"Color-coding: AGE\" dialog")));
      await session.step(52, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(53, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(54, "When user selects \"Off\" in Type input in \"Color-coding: AGE\" dialog", () => selectIn(page, "Off", el("Type input in \"Color-coding: AGE\" dialog")));
      await session.step(55, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(56, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(57, "When user selects \"Linear\" in Type input in \"Color-coding: AGE\" dialog", () => selectIn(page, "Linear", el("Type input in \"Color-coding: AGE\" dialog")));
      await session.step(58, "Then \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
      await session.step(59, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(60, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(61, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(62, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(63, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Apply to text sets the column's text flag", async () => {
      await session.step(67, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(68, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(69, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(70, "When user selects \"text\" in Apply to input in \"Color-coding: AGE\" dialog", () => selectIn(page, "text", el("Apply to input in \"Color-coding: AGE\" dialog")));
      await session.step(71, "Then the text of \"AGE\" column should be color-coded", () => textColorCoded(page, "AGE"));
      await session.step(72, "When user selects \"background\" in Apply to input in \"Color-coding: AGE\" dialog", () => selectIn(page, "background", el("Apply to input in \"Color-coding: AGE\" dialog")));
      await session.step(73, "And user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(74, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(75, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(76, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Min and Max clamp the scheme", async () => {
      await session.step(80, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(81, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(82, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(83, "And the \"color of cell 1 of AGE\" and \"color of cell 2 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 2 of AGE", el("grid")));
      await session.step(84, "When user enters \"50\" into Min input in \"Color-coding: AGE\" dialog", () => enterInto(page, "50", el("Min input in \"Color-coding: AGE\" dialog")));
      await session.step(85, "And user enters \"60\" into Max input in \"Color-coding: AGE\" dialog", () => enterInto(page, "60", el("Max input in \"Color-coding: AGE\" dialog")));
      await session.step(86, "Then the \"color of cell 1 of AGE\" and \"color of cell 2 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 2 of AGE", el("grid")));
      await session.step(87, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(88, "When user clears Min input in \"Color-coding: AGE\" dialog", () => clearField(page, el("Min input in \"Color-coding: AGE\" dialog")));
      await session.step(89, "And user clears Max input in \"Color-coding: AGE\" dialog", () => clearField(page, el("Max input in \"Color-coding: AGE\" dialog")));
      await session.step(90, "Then the \"color of cell 1 of AGE\" and \"color of cell 2 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 2 of AGE", el("grid")));
      await session.step(91, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(92, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(93, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(94, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The invert icon swaps the ends of the gradient", async () => {
      await session.step(98, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(99, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(100, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(101, "When user remembers the \"color of cell 3 of AGE\" reading of grid", () => rememberReading(page, "color of cell 3 of AGE", el("grid")));
      await session.step(102, "And user clicks on \"exchange\" icon in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"exchange\" icon in \"Color-coding: AGE\" dialog")));
      await session.step(103, "Then the color scheme of \"AGE\" column should be \"#FF0000, #FFFFFF, #0000FF\"", () => colorSchemeIs(page, "AGE", "#FF0000, #FFFFFF, #0000FF"));
      await session.step(104, "And the \"color of cell 3 of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(105, "When user clicks on \"exchange\" icon in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"exchange\" icon in \"Color-coding: AGE\" dialog")));
      await session.step(106, "Then the color scheme of \"AGE\" column should be \"#0000FF, #FFFFFF, #FF0000\"", () => colorSchemeIs(page, "AGE", "#0000FF, #FFFFFF, #FF0000"));
      await session.step(107, "And the \"color of cell 3 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(108, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(109, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(110, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(111, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit scheme: absolute values ignore Min and Max, a stop is removed, and the edited scheme survives an invert", async () => {
      await session.step(115, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(116, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(117, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(118, "And the color scheme of \"AGE\" column should be \"#0000FF, #FFFFFF, #FF0000\"", () => colorSchemeIs(page, "AGE", "#0000FF, #FFFFFF, #FF0000"));
      await session.step(119, "When user enters \"50\" into Min input in \"Color-coding: AGE\" dialog", () => enterInto(page, "50", el("Min input in \"Color-coding: AGE\" dialog")));
      await session.step(120, "And user enters \"60\" into Max input in \"Color-coding: AGE\" dialog", () => enterInto(page, "60", el("Max input in \"Color-coding: AGE\" dialog")));
      await session.step(121, "And user clicks on \"Edit scheme...\" action in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Edit scheme...\" action in \"Color-coding: AGE\" dialog")));
      await session.step(122, "Then \"Edit color scheme\" dialog should be visible", () => shouldBe(page, el("\"Edit color scheme\" dialog"), "visible"));
      await session.step(123, "When user remembers the \"color of cell 3 of AGE\" reading of grid", () => rememberReading(page, "color of cell 3 of AGE", el("grid")));
      await session.step(124, "And user checks \"Use absolute values\" input in \"Edit color scheme\" dialog", () => check(page, el("\"Use absolute values\" input in \"Edit color scheme\" dialog")));
      await session.step(125, "Then \"AGE\" column should have tag \".%color-coding-linear-is-absolute\" equal to \"true\"", () => columnTag(page, "AGE", ".%color-coding-linear-is-absolute", "true"));
      await session.step(126, "And the \"color of cell 3 of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(127, "When user unchecks \"Use absolute values\" input in \"Edit color scheme\" dialog", () => uncheck(page, el("\"Use absolute values\" input in \"Edit color scheme\" dialog")));
      await session.step(128, "Then the \"color of cell 3 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(129, "When user remembers the \"color of cell 3 of AGE\" reading of grid", () => rememberReading(page, "color of cell 3 of AGE", el("grid")));
      await session.step(130, "And user clicks on first Remove row button in \"Edit color scheme\" dialog", () => clickOn(page, el("first Remove row button in \"Edit color scheme\" dialog")));
      await session.step(131, "Then the color scheme of \"AGE\" column should be \"#FFFFFF, #FF0000\"", () => colorSchemeIs(page, "AGE", "#FFFFFF, #FF0000"));
      await session.step(132, "And the \"color of cell 3 of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(133, "When user clicks on OK button in \"Edit color scheme\" dialog", () => clickOn(page, el("OK button in \"Edit color scheme\" dialog")));
      await session.step(134, "Then \"Edit color scheme\" dialog should be hidden", () => shouldBe(page, el("\"Edit color scheme\" dialog"), "hidden"));
      await session.step(135, "When user clicks on \"exchange\" icon in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"exchange\" icon in \"Color-coding: AGE\" dialog")));
      await session.step(136, "Then the color scheme of \"AGE\" column should be \"#FF0000, #FFFFFF\"", () => colorSchemeIs(page, "AGE", "#FF0000, #FFFFFF"));
      await session.step(137, "When user clears Min input in \"Color-coding: AGE\" dialog", () => clearField(page, el("Min input in \"Color-coding: AGE\" dialog")));
      await session.step(138, "And user clears Max input in \"Color-coding: AGE\" dialog", () => clearField(page, el("Max input in \"Color-coding: AGE\" dialog")));
      await session.step(139, "And user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(140, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(141, "And \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
      await session.step(142, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(143, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(144, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
