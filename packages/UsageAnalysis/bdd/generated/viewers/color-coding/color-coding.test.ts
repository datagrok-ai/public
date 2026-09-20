/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/color-coding/color-coding.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, categoricalColorIs, colorAgain, colorCategorical, colorCodedAs, colorCodedCategorically, colorConditional, colorInverted, colorLinear, colorLinearThrough, colorLinked, colorLinkedText, colorLinkedTo, colorOff, colorPickUp, colorSchemeIs, noColorCoding, removeColumn, textColorCoded} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaColor, loadLayout, noErrors, readingDoesNotRead, readingNotAsRemembered, readingReads, readingsDiffer, readingsEqual, rememberReading, repainted, saveLayoutToServer, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Colour coding across columns", () => {
  const session = feature(test, "features/viewers/color-coding/color-coding.feature", import.meta.url);
  test("Colour coding across columns", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 12, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(24, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(25, "And \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
    await session.step(26, "And \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await run.scenario("A colouring switched off and back on keeps the colours it was given", async () => {
      await session.step(29, "When user colors \"SEX\" column categorically:", () => colorCategorical(page, "SEX", [["M","#3366CC"],["F","#CC6699"]]));
      await session.step(32, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(33, "And the categorical color of \"M\" in \"SEX\" column should be \"3366CC\"", () => categoricalColorIs(page, "M", "SEX", "3366CC"));
      await session.step(34, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(35, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(36, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(37, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(38, "When user colors \"SEX\" column categorically again", () => colorAgain(page, "SEX", "categorically"));
      await session.step(39, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(40, "And the categorical color of \"M\" in \"SEX\" column should be \"3366CC\"", () => categoricalColorIs(page, "M", "SEX", "3366CC"));
      await session.step(41, "And the categorical color of \"F\" in \"SEX\" column should be \"CC6699\"", () => categoricalColorIs(page, "F", "SEX", "CC6699"));
      await session.step(42, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(43, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(44, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Applying one column's colouring to another gives it the same colours", async () => {
      await session.step(47, "When user colors \"AGE\" column conditionally:", () => colorConditional(page, "AGE", [["<30","#00FF00"],["30-60","#FFFF00"],[">60","#FF0000"]]));
      await session.step(51, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(52, "And the \"color of cell 1 of AGE\" reading of grid should be \"#00ff00\"", () => readingReads(page, "color of cell 1 of AGE", el("grid"), "#00ff00"));
      await session.step(53, "And the \"color of cell 3 of AGE\" reading of grid should be \"#ffff00\"", () => readingReads(page, "color of cell 3 of AGE", el("grid"), "#ffff00"));
      await session.step(54, "When user adds a calculated column \"Age_copy\" with formula \"${AGE}\"", () => addCalculated(page, "Age_copy", "${AGE}"));
      await session.step(55, "Then \"Age_copy\" column should have no color coding", () => noColorCoding(page, "Age_copy"));
      await session.step(56, "And the \"color of cell 1 of AGE\" and \"color of cell 1 of Age_copy\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 1 of Age_copy", el("grid")));
      await session.step(57, "When user applies the coloring of \"AGE\" column to \"Age_copy\" column", () => colorPickUp(page, "AGE", "Age_copy"));
      await session.step(58, "Then \"Age_copy\" column should be color-coded conditionally", () => colorCodedAs(page, "Age_copy", "conditionally"));
      await session.step(59, "And the \"color of cell 1 of AGE\" and \"color of cell 1 of Age_copy\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 1 of Age_copy", el("grid")));
      await session.step(60, "And the \"color of cell 3 of AGE\" and \"color of cell 3 of Age_copy\" readings of grid should be the same", () => readingsEqual(page, "color of cell 3 of AGE", "color of cell 3 of Age_copy", el("grid")));
      await session.step(61, "When user removes \"Age_copy\" column", () => removeColumn(page, "Age_copy"));
      await session.step(62, "Then the table should not have a column \"Age_copy\"", () => hasNoColumn(page, "Age_copy"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A linked column paints its cells with the source column's colours", async () => {
      await session.step(66, "When user colors \"RACE\" column linked to \"AGE\" column", () => colorLinked(page, "RACE", "AGE"));
      await session.step(67, "Then \"RACE\" column should be color-coded linked", () => colorCodedAs(page, "RACE", "linked"));
      await session.step(68, "And the coloring of \"RACE\" column should be linked to \"AGE\" column", () => colorLinkedTo(page, "RACE", "AGE"));
      await session.step(69, "And the \"color of cell 1 of RACE\" reading of grid should be \"#00ff00\"", () => readingReads(page, "color of cell 1 of RACE", el("grid"), "#00ff00"));
      await session.step(70, "And the \"color of cell 3 of RACE\" reading of grid should be \"#ffff00\"", () => readingReads(page, "color of cell 3 of RACE", el("grid"), "#ffff00"));
      await session.step(71, "And the \"color of cell 1 of RACE\" and \"color of cell 3 of RACE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of RACE", "color of cell 3 of RACE", el("grid")));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A change of the source's type leaves the link in place", async () => {
      await session.step(75, "When user colors \"AGE\" column linearly from \"#0000FF\" to \"#FF0000\"", () => colorLinear(page, "AGE", "#0000FF", "#FF0000"));
      await session.step(76, "Then \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
      await session.step(77, "And \"RACE\" column should be color-coded linked", () => colorCodedAs(page, "RACE", "linked"));
      await session.step(78, "And the coloring of \"RACE\" column should be linked to \"AGE\" column", () => colorLinkedTo(page, "RACE", "AGE"));
      await session.step(79, "And the \"color of cell 1 of RACE\" and \"color of cell 1 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of RACE", "color of cell 1 of AGE", el("grid")));
      await session.step(80, "And the \"color of cell 1 of RACE\" reading of grid should not be \"#00ff00\"", () => readingDoesNotRead(page, "color of cell 1 of RACE", el("grid"), "#00ff00"));
      await session.step(81, "When user colors \"AGE\" column conditionally:", () => colorConditional(page, "AGE", [["<40","#00FF00"],[">40","#FF0000"]]));
      await session.step(84, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(85, "And \"RACE\" column should be color-coded linked", () => colorCodedAs(page, "RACE", "linked"));
      await session.step(86, "And the \"color of cell 1 of RACE\" reading of grid should be \"#00ff00\"", () => readingReads(page, "color of cell 1 of RACE", el("grid"), "#00ff00"));
      await session.step(87, "And the \"color of cell 3 of RACE\" reading of grid should be \"#ff0000\"", () => readingReads(page, "color of cell 3 of RACE", el("grid"), "#ff0000"));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A linked colouring applied to the text paints the letters, not the cell", async () => {
      await session.step(91, "When user colors the text of \"HEIGHT\" column linked to \"AGE\" column", () => colorLinkedText(page, "HEIGHT", "AGE"));
      await session.step(92, "Then \"HEIGHT\" column should be color-coded linked", () => colorCodedAs(page, "HEIGHT", "linked"));
      await session.step(93, "And the text of \"HEIGHT\" column should be color-coded", () => textColorCoded(page, "HEIGHT"));
      await session.step(94, "And the \"color of cell 1 of HEIGHT\" and \"color of cell 1 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of HEIGHT", "color of cell 1 of AGE", el("grid")));
      await session.step(95, "And the \"color of cell 1 of HEIGHT\" and \"color of cell 3 of HEIGHT\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of HEIGHT", "color of cell 3 of HEIGHT", el("grid")));
      await session.step(96, "And the \"cell 1 of HEIGHT\" area of grid should contain the color \"#00FF00\"", () => areaColor(page, "cell 1 of HEIGHT", el("grid"), "#00FF00"));
      await session.step(97, "When user removes the coloring of \"HEIGHT\" column", () => colorOff(page, "HEIGHT"));
      await session.step(98, "And user removes the coloring of \"RACE\" column", () => colorOff(page, "RACE"));
      await session.step(99, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A five-level chain of links all reports Linked", async () => {
      await session.step(103, "When user colors \"SEX\" column linked to \"AGE\" column", () => colorLinked(page, "SEX", "AGE"));
      await session.step(104, "And user colors \"DIS_POP\" column linked to \"SEX\" column", () => colorLinked(page, "DIS_POP", "SEX"));
      await session.step(105, "And user colors \"CONTROL\" column linked to \"DIS_POP\" column", () => colorLinked(page, "CONTROL", "DIS_POP"));
      await session.step(106, "And user colors \"STARTED\" column linked to \"CONTROL\" column", () => colorLinked(page, "STARTED", "CONTROL"));
      await session.step(107, "Then \"SEX\" column should be color-coded linked", () => colorCodedAs(page, "SEX", "linked"));
      await session.step(108, "And \"DIS_POP\" column should be color-coded linked", () => colorCodedAs(page, "DIS_POP", "linked"));
      await session.step(109, "And \"CONTROL\" column should be color-coded linked", () => colorCodedAs(page, "CONTROL", "linked"));
      await session.step(110, "And \"STARTED\" column should be color-coded linked", () => colorCodedAs(page, "STARTED", "linked"));
      await session.step(111, "And the coloring of \"STARTED\" column should be linked to \"CONTROL\" column", () => colorLinkedTo(page, "STARTED", "CONTROL"));
      await session.step(112, "And the \"color of cell 1 of SEX\" reading of grid should be \"#00ff00\"", () => readingReads(page, "color of cell 1 of SEX", el("grid"), "#00ff00"));
      await session.step(113, "And the \"color of cell 3 of SEX\" reading of grid should be \"#ff0000\"", () => readingReads(page, "color of cell 3 of SEX", el("grid"), "#ff0000"));
      await session.step(114, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(115, "And user removes the coloring of \"DIS_POP\" column", () => colorOff(page, "DIS_POP"));
      await session.step(116, "And user removes the coloring of \"CONTROL\" column", () => colorOff(page, "CONTROL"));
      await session.step(117, "And user removes the coloring of \"STARTED\" column", () => colorOff(page, "STARTED"));
      await session.step(118, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Inverting a linear scheme swaps the ends of the gradient", async () => {
      await session.step(121, "When user colors \"AGE\" column linearly from \"#0000FF\" to \"#FF0000\"", () => colorLinear(page, "AGE", "#0000FF", "#FF0000"));
      await session.step(122, "Then the color scheme of \"AGE\" column should be \"#0000FF, #FF0000\"", () => colorSchemeIs(page, "AGE", "#0000FF, #FF0000"));
      await session.step(123, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(124, "When user remembers the \"color of cell 1 of AGE\" reading of grid", () => rememberReading(page, "color of cell 1 of AGE", el("grid")));
      await session.step(125, "And user inverts the color scheme of \"AGE\" column", () => colorInverted(page, "AGE"));
      await session.step(126, "Then the color scheme of \"AGE\" column should be \"#FF0000, #0000FF\"", () => colorSchemeIs(page, "AGE", "#FF0000, #0000FF"));
      await session.step(127, "And grid should have repainted", () => repainted(page, el("grid")));
      await session.step(128, "And the \"color of cell 1 of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "color of cell 1 of AGE", el("grid")));
      await session.step(129, "When user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
      await session.step(130, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(131, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A boolean column takes the default categorical colours, a date column a three-stop scheme", async () => {
      await session.step(134, "When user colors \"CONTROL\" column categorically again", () => colorAgain(page, "CONTROL", "categorically"));
      await session.step(135, "Then \"CONTROL\" column should be color-coded categorically", () => colorCodedCategorically(page, "CONTROL"));
      await session.step(136, "And the \"color of cell 1 of CONTROL\" and \"color of cell 3 of CONTROL\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of CONTROL", "color of cell 3 of CONTROL", el("grid")));
      await session.step(137, "When user colors \"STARTED\" column linearly through \"#0000FF, #FFFFFF, #FF0000\"", () => colorLinearThrough(page, "STARTED", "#0000FF, #FFFFFF, #FF0000"));
      await session.step(138, "Then \"STARTED\" column should be color-coded linearly", () => colorCodedAs(page, "STARTED", "linearly"));
      await session.step(139, "And the color scheme of \"STARTED\" column should be \"#0000FF, #FFFFFF, #FF0000\"", () => colorSchemeIs(page, "STARTED", "#0000FF, #FFFFFF, #FF0000"));
      await session.step(140, "And the \"color of cell 1 of STARTED\" and \"color of cell 2 of STARTED\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of STARTED", "color of cell 2 of STARTED", el("grid")));
      await session.step(141, "When user removes the coloring of \"CONTROL\" column", () => colorOff(page, "CONTROL"));
      await session.step(142, "Then \"CONTROL\" column should have no color coding", () => noColorCoding(page, "CONTROL"));
      await session.step(143, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("AGE and STARTED switched off and back on keep the schemes they were given", async () => {
      await session.step(146, "When user colors \"AGE\" column conditionally:", () => colorConditional(page, "AGE", [["<30","#00FF00"],["30-60","#FFFF00"],[">60","#FF0000"]]));
      await session.step(150, "Then the \"color of cell 1 of AGE\" reading of grid should be \"#00ff00\"", () => readingReads(page, "color of cell 1 of AGE", el("grid"), "#00ff00"));
      await session.step(151, "When user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
      await session.step(152, "And user removes the coloring of \"STARTED\" column", () => colorOff(page, "STARTED"));
      await session.step(153, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(154, "And \"STARTED\" column should have no color coding", () => noColorCoding(page, "STARTED"));
      await session.step(155, "And the \"color of cell 1 of STARTED\" and \"color of cell 2 of STARTED\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of STARTED", "color of cell 2 of STARTED", el("grid")));
      await session.step(156, "When user colors \"AGE\" column conditionally again", () => colorAgain(page, "AGE", "conditionally"));
      await session.step(157, "And user colors \"STARTED\" column linearly again", () => colorAgain(page, "STARTED", "linearly"));
      await session.step(158, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(159, "And the \"color of cell 1 of AGE\" reading of grid should be \"#00ff00\"", () => readingReads(page, "color of cell 1 of AGE", el("grid"), "#00ff00"));
      await session.step(160, "And the \"color of cell 3 of AGE\" reading of grid should be \"#ffff00\"", () => readingReads(page, "color of cell 3 of AGE", el("grid"), "#ffff00"));
      await session.step(161, "And \"STARTED\" column should be color-coded linearly", () => colorCodedAs(page, "STARTED", "linearly"));
      await session.step(162, "And the color scheme of \"STARTED\" column should be \"#0000FF, #FFFFFF, #FF0000\"", () => colorSchemeIs(page, "STARTED", "#0000FF, #FFFFFF, #FF0000"));
      await session.step(163, "And the \"color of cell 1 of STARTED\" and \"color of cell 2 of STARTED\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of STARTED", "color of cell 2 of STARTED", el("grid")));
      await session.step(164, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pick Up / Apply copies a categorical colouring and a linear scheme to other columns", async () => {
      await session.step(167, "When user adds a calculated column \"Race_copy\" with formula \"${RACE}\"", () => addCalculated(page, "Race_copy", "${RACE}"));
      await session.step(168, "Then \"Race_copy\" column should have no color coding", () => noColorCoding(page, "Race_copy"));
      await session.step(169, "And \"HEIGHT\" column should have no color coding", () => noColorCoding(page, "HEIGHT"));
      await session.step(170, "When user colors \"RACE\" column categorically:", () => colorCategorical(page, "RACE", [["Caucasian","#3366CC"],["Other","#CC6699"]]));
      await session.step(173, "And user applies the coloring of \"RACE\" column to \"Race_copy\" column", () => colorPickUp(page, "RACE", "Race_copy"));
      await session.step(174, "Then \"Race_copy\" column should be color-coded categorically", () => colorCodedCategorically(page, "Race_copy"));
      await session.step(175, "And the categorical color of \"Caucasian\" in \"Race_copy\" column should be \"3366CC\"", () => categoricalColorIs(page, "Caucasian", "Race_copy", "3366CC"));
      await session.step(176, "And the categorical color of \"Other\" in \"Race_copy\" column should be \"CC6699\"", () => categoricalColorIs(page, "Other", "Race_copy", "CC6699"));
      await session.step(177, "And the \"color of cell 1 of Race_copy\" and \"color of cell 1 of RACE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of Race_copy", "color of cell 1 of RACE", el("grid")));
      await session.step(178, "When user applies the coloring of \"STARTED\" column to \"HEIGHT\" column", () => colorPickUp(page, "STARTED", "HEIGHT"));
      await session.step(179, "Then \"HEIGHT\" column should be color-coded linearly", () => colorCodedAs(page, "HEIGHT", "linearly"));
      await session.step(180, "And the color scheme of \"HEIGHT\" column should be \"#0000FF, #FFFFFF, #FF0000\"", () => colorSchemeIs(page, "HEIGHT", "#0000FF, #FFFFFF, #FF0000"));
      await session.step(181, "And the \"color of cell 1 of HEIGHT\" and \"color of cell 2 of HEIGHT\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of HEIGHT", "color of cell 2 of HEIGHT", el("grid")));
      await session.step(182, "When user removes the coloring of \"HEIGHT\" column", () => colorOff(page, "HEIGHT"));
      await session.step(183, "And user removes the coloring of \"RACE\" column", () => colorOff(page, "RACE"));
      await session.step(184, "And user removes the coloring of \"STARTED\" column", () => colorOff(page, "STARTED"));
      await session.step(185, "And user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
      await session.step(186, "And user removes \"Race_copy\" column", () => removeColumn(page, "Race_copy"));
      await session.step(187, "Then the table should not have a column \"Race_copy\"", () => hasNoColumn(page, "Race_copy"));
      await session.step(188, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column's colouring comes back from a saved layout on the reopened table", async () => {
      await session.step(191, "When user colors \"SEX\" column categorically:", () => colorCategorical(page, "SEX", [["M","#3366CC"],["F","#CC6699"]]));
      await session.step(194, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(195, "And user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(196, "And user closes all views", () => closeAllViews(page));
      await session.step(197, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(198, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(199, "When user loads the saved layout", () => loadLayout(page));
      await session.step(200, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(201, "And \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(202, "And the categorical color of \"M\" in \"SEX\" column should be \"3366CC\"", () => categoricalColorIs(page, "M", "SEX", "3366CC"));
      await session.step(203, "And the categorical color of \"F\" in \"SEX\" column should be \"CC6699\"", () => categoricalColorIs(page, "F", "SEX", "CC6699"));
      await session.step(204, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(205, "When user clicks on close icon of scatter plot viewer", () => clickOn(page, el("close icon of scatter plot viewer")));
      await session.step(206, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(207, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(208, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Linked colourings come back from a saved layout and from a saved project", async () => {
      await session.step(211, "When user colors \"AGE\" column conditionally:", () => colorConditional(page, "AGE", [["<40","#00FF00"],[">40","#FF0000"]]));
      await session.step(214, "And user colors \"RACE\" column linked to \"AGE\" column", () => colorLinked(page, "RACE", "AGE"));
      await session.step(215, "And user colors \"HEIGHT\" column linked to \"RACE\" column", () => colorLinked(page, "HEIGHT", "RACE"));
      await session.step(216, "Then the coloring of \"RACE\" column should be linked to \"AGE\" column", () => colorLinkedTo(page, "RACE", "AGE"));
      await session.step(217, "And the coloring of \"HEIGHT\" column should be linked to \"RACE\" column", () => colorLinkedTo(page, "HEIGHT", "RACE"));
      await session.step(218, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(219, "And user removes the coloring of \"HEIGHT\" column", () => colorOff(page, "HEIGHT"));
      await session.step(220, "And user removes the coloring of \"RACE\" column", () => colorOff(page, "RACE"));
      await session.step(221, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
      await session.step(222, "And \"HEIGHT\" column should have no color coding", () => noColorCoding(page, "HEIGHT"));
      await session.step(223, "When user loads the saved layout", () => loadLayout(page));
      await session.step(224, "Then \"RACE\" column should be color-coded linked", () => colorCodedAs(page, "RACE", "linked"));
      await session.step(225, "And the coloring of \"RACE\" column should be linked to \"AGE\" column", () => colorLinkedTo(page, "RACE", "AGE"));
      await session.step(226, "And \"HEIGHT\" column should be color-coded linked", () => colorCodedAs(page, "HEIGHT", "linked"));
      await session.step(227, "And the coloring of \"HEIGHT\" column should be linked to \"RACE\" column", () => colorLinkedTo(page, "HEIGHT", "RACE"));
      await session.step(228, "And the \"color of cell 1 of RACE\" reading of grid should be \"#00ff00\"", () => readingReads(page, "color of cell 1 of RACE", el("grid"), "#00ff00"));
      await session.step(229, "And the \"color of cell 3 of RACE\" reading of grid should be \"#ff0000\"", () => readingReads(page, "color of cell 3 of RACE", el("grid"), "#ff0000"));
      await session.step(230, "When user saves the current view as project \"bdd color coding links\"", () => saveAsProject(page, "bdd color coding links"));
      await session.step(231, "And user closes all views", () => closeAllViews(page));
      await session.step(232, "And user opens the \"bdd color coding links\" project", () => openProject(page, "bdd color coding links"));
      await session.step(233, "Then \"RACE\" column should be color-coded linked", () => colorCodedAs(page, "RACE", "linked"));
      await session.step(234, "And the coloring of \"RACE\" column should be linked to \"AGE\" column", () => colorLinkedTo(page, "RACE", "AGE"));
      await session.step(235, "And \"HEIGHT\" column should be color-coded linked", () => colorCodedAs(page, "HEIGHT", "linked"));
      await session.step(236, "And the coloring of \"HEIGHT\" column should be linked to \"RACE\" column", () => colorLinkedTo(page, "HEIGHT", "RACE"));
      await session.step(237, "And the \"color of cell 1 of RACE\" reading of grid should be \"#00ff00\"", () => readingReads(page, "color of cell 1 of RACE", el("grid"), "#00ff00"));
      await session.step(238, "And the \"color of cell 3 of RACE\" reading of grid should be \"#ff0000\"", () => readingReads(page, "color of cell 3 of RACE", el("grid"), "#ff0000"));
      await session.step(239, "When user removes the coloring of \"HEIGHT\" column", () => colorOff(page, "HEIGHT"));
      await session.step(240, "And user removes the coloring of \"RACE\" column", () => colorOff(page, "RACE"));
      await session.step(241, "And user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
      await session.step(242, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
      await session.step(243, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
