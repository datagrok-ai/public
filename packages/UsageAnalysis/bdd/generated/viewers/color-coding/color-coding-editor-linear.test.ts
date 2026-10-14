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
import {check, clickOn, enterInto, selectIn, shouldBe, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorCodedAs, colorSchemeIs, noColorCoding, textColorCoded} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromAreaContextMenu, readingAsRemembered, readingNotAsRemembered, readingsDiffer, readingsEqual, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
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
      await session.step(27, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
      await session.step(28, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(29, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(30, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(31, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(32, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Type choice of the editor switches the column's coding", async () => {
      await session.step(36, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(37, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")));
      await session.step(38, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(39, "And Type input in \"Color-coding: AGE\" dialog should have value \"Linear\"", () => shouldHaveValue(page, el("Type input in \"Color-coding: AGE\" dialog"), "Linear"));
      await session.step(40, "When user selects \"Conditional\" in Type input in \"Color-coding: AGE\" dialog", () => selectIn(page, "Conditional", el("Type input in \"Color-coding: AGE\" dialog")));
      await session.step(41, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(42, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(43, "When user selects \"Off\" in Type input in \"Color-coding: AGE\" dialog", () => selectIn(page, "Off", el("Type input in \"Color-coding: AGE\" dialog")));
      await session.step(44, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(45, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(46, "When user selects \"Linear\" in Type input in \"Color-coding: AGE\" dialog", () => selectIn(page, "Linear", el("Type input in \"Color-coding: AGE\" dialog")));
      await session.step(47, "Then \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
      await session.step(48, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Apply to text sets the column's text flag", async () => {
      await session.step(52, "When user selects \"text\" in Apply to input in \"Color-coding: AGE\" dialog", () => selectIn(page, "text", el("Apply to input in \"Color-coding: AGE\" dialog")));
      await session.step(53, "Then the text of \"AGE\" column should be color-coded", () => textColorCoded(page, "AGE"));
      await session.step(54, "When user selects \"background\" in Apply to input in \"Color-coding: AGE\" dialog", () => selectIn(page, "background", el("Apply to input in \"Color-coding: AGE\" dialog")));
      await session.step(55, "Then the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Min and Max clamp the scheme", async () => {
      await session.step(59, "Then the \"color of cell 1 of AGE\" and \"color of cell 2 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 2 of AGE", el("grid")));
      await session.step(60, "When user enters \"50\" into Min input in \"Color-coding: AGE\" dialog", () => enterInto(page, "50", el("Min input in \"Color-coding: AGE\" dialog")));
      await session.step(61, "And user enters \"60\" into Max input in \"Color-coding: AGE\" dialog", () => enterInto(page, "60", el("Max input in \"Color-coding: AGE\" dialog")));
      await session.step(62, "Then \"AGE\" column should have tag \".color-coding-scheme-min\" equal to \"50\"", () => columnTag(page, "AGE", ".color-coding-scheme-min", "50"));
      await session.step(63, "And \"AGE\" column should have tag \".color-coding-scheme-max\" equal to \"60\"", () => columnTag(page, "AGE", ".color-coding-scheme-max", "60"));
      await session.step(64, "And the \"color of cell 1 of AGE\" and \"color of cell 2 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 2 of AGE", el("grid")));
      await session.step(65, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The invert icon swaps the ends of the gradient", async () => {
      await session.step(69, "When user remembers the \"color of cell 3 of AGE\" reading of grid", () => rememberReading(page, "color of cell 3 of AGE", el("grid")));
      await session.step(70, "And user clicks on \"exchange\" icon in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"exchange\" icon in \"Color-coding: AGE\" dialog")));
      await session.step(71, "Then the color scheme of \"AGE\" column should be \"#FF0000, #FFFFFF, #0000FF\"", () => colorSchemeIs(page, "AGE", "#FF0000, #FFFFFF, #0000FF"));
      await session.step(72, "And the \"color of cell 3 of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(73, "When user clicks on \"exchange\" icon in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"exchange\" icon in \"Color-coding: AGE\" dialog")));
      await session.step(74, "Then the color scheme of \"AGE\" column should be \"#0000FF, #FFFFFF, #FF0000\"", () => colorSchemeIs(page, "AGE", "#0000FF, #FFFFFF, #FF0000"));
      await session.step(75, "And the \"color of cell 3 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit scheme removes a stop, and the edited scheme survives an invert", async () => {
      await session.step(79, "Then the color scheme of \"AGE\" column should be \"#0000FF, #FFFFFF, #FF0000\"", () => colorSchemeIs(page, "AGE", "#0000FF, #FFFFFF, #FF0000"));
      await session.step(80, "When user clicks on \"Edit scheme...\" action in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Edit scheme...\" action in \"Color-coding: AGE\" dialog")));
      await session.step(81, "Then \"Edit color scheme\" dialog should be visible", () => shouldBe(page, el("\"Edit color scheme\" dialog"), "visible"));
      await session.step(82, "When user remembers the \"color of cell 3 of AGE\" reading of grid", () => rememberReading(page, "color of cell 3 of AGE", el("grid")));
      await session.step(83, "And user checks \"Use absolute values\" input in \"Edit color scheme\" dialog", () => check(page, el("\"Use absolute values\" input in \"Edit color scheme\" dialog")));
      await session.step(84, "Then \"AGE\" column should have tag \".%color-coding-linear-is-absolute\" equal to \"true\"", () => columnTag(page, "AGE", ".%color-coding-linear-is-absolute", "true"));
      await session.step(85, "And the \"color of cell 3 of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(86, "When user unchecks \"Use absolute values\" input in \"Edit color scheme\" dialog", () => uncheck(page, el("\"Use absolute values\" input in \"Edit color scheme\" dialog")));
      await session.step(87, "Then the \"color of cell 3 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(88, "When user remembers the \"color of cell 3 of AGE\" reading of grid", () => rememberReading(page, "color of cell 3 of AGE", el("grid")));
      await session.step(89, "And user clicks on first Remove row button in \"Edit color scheme\" dialog", () => clickOn(page, el("first Remove row button in \"Edit color scheme\" dialog")));
      await session.step(90, "Then the color scheme of \"AGE\" column should be \"#FFFFFF, #FF0000\"", () => colorSchemeIs(page, "AGE", "#FFFFFF, #FF0000"));
      await session.step(91, "And the \"color of cell 3 of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "color of cell 3 of AGE", el("grid")));
      await session.step(92, "When user clicks on OK button in \"Edit color scheme\" dialog", () => clickOn(page, el("OK button in \"Edit color scheme\" dialog")));
      await session.step(93, "Then \"Edit color scheme\" dialog should be hidden", () => shouldBe(page, el("\"Edit color scheme\" dialog"), "hidden"));
      await session.step(94, "When user clicks on \"exchange\" icon in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"exchange\" icon in \"Color-coding: AGE\" dialog")));
      await session.step(95, "Then the color scheme of \"AGE\" column should be \"#FF0000, #FFFFFF\"", () => colorSchemeIs(page, "AGE", "#FF0000, #FFFFFF"));
      await session.step(96, "When user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
      await session.step(97, "Then \"Color-coding: AGE\" dialog should be hidden", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "hidden"));
      await session.step(98, "And \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
      await session.step(99, "When user picks \"Color Coding > Off\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Off", "header AGE", el("grid")));
      await session.step(100, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
