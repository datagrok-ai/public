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
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, categoricalColorIs, colorAgain, colorCategorical, colorCodedAs, colorCodedCategorically, colorConditional, colorInverted, colorLinear, colorLinked, colorLinkedText, colorLinkedTo, colorOff, colorPickUp, colorSchemeIs, noColorCoding, removeColumn, textColorCoded} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColor, noErrors, readingDoesNotRead, readingNotAsRemembered, readingReads, readingsDiffer, readingsEqual, rememberReading, repainted, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Colour coding across columns", () => {
  const session = feature(test, "features/viewers/color-coding/color-coding.feature", import.meta.url);
  test("Colour coding across columns", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(16, "And \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
    await session.step(17, "And \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
    await run.scenario("A colouring switched off and back on keeps the colours it was given", async () => {
      await session.step(20, "When user colors \"SEX\" column categorically:", () => colorCategorical(page, "SEX", [["M","#3366CC"],["F","#CC6699"]]));
      await session.step(23, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(24, "And the categorical color of \"M\" in \"SEX\" column should be \"3366CC\"", () => categoricalColorIs(page, "M", "SEX", "3366CC"));
      await session.step(25, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(26, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(27, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(28, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(29, "When user colors \"SEX\" column categorically again", () => colorAgain(page, "SEX", "categorically"));
      await session.step(30, "Then \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(31, "And the categorical color of \"M\" in \"SEX\" column should be \"3366CC\"", () => categoricalColorIs(page, "M", "SEX", "3366CC"));
      await session.step(32, "And the categorical color of \"F\" in \"SEX\" column should be \"CC6699\"", () => categoricalColorIs(page, "F", "SEX", "CC6699"));
      await session.step(33, "And the \"color of cell 1 of SEX\" and \"color of cell 4 of SEX\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of SEX", "color of cell 4 of SEX", el("grid")));
      await session.step(34, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(35, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Applying one column's colouring to another gives it the same colours", async () => {
      await session.step(38, "When user colors \"AGE\" column conditionally:", () => colorConditional(page, "AGE", [["<30","#00FF00"],["30-60","#FFFF00"],[">60","#FF0000"]]));
      await session.step(42, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(43, "And the \"color of cell 1 of AGE\" reading of grid should be \"#00ff00\"", () => readingReads(page, "color of cell 1 of AGE", el("grid"), "#00ff00"));
      await session.step(44, "And the \"color of cell 3 of AGE\" reading of grid should be \"#ffff00\"", () => readingReads(page, "color of cell 3 of AGE", el("grid"), "#ffff00"));
      await session.step(45, "When user adds a calculated column \"Age_copy\" with formula \"${AGE}\"", () => addCalculated(page, "Age_copy", "${AGE}"));
      await session.step(46, "Then \"Age_copy\" column should have no color coding", () => noColorCoding(page, "Age_copy"));
      await session.step(47, "And the \"color of cell 1 of AGE\" and \"color of cell 1 of Age_copy\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 1 of Age_copy", el("grid")));
      await session.step(48, "When user applies the coloring of \"AGE\" column to \"Age_copy\" column", () => colorPickUp(page, "AGE", "Age_copy"));
      await session.step(49, "Then \"Age_copy\" column should be color-coded conditionally", () => colorCodedAs(page, "Age_copy", "conditionally"));
      await session.step(50, "And the \"color of cell 1 of AGE\" and \"color of cell 1 of Age_copy\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of AGE", "color of cell 1 of Age_copy", el("grid")));
      await session.step(51, "And the \"color of cell 3 of AGE\" and \"color of cell 3 of Age_copy\" readings of grid should be the same", () => readingsEqual(page, "color of cell 3 of AGE", "color of cell 3 of Age_copy", el("grid")));
      await session.step(52, "When user removes \"Age_copy\" column", () => removeColumn(page, "Age_copy"));
      await session.step(53, "Then the table should not have a column \"Age_copy\"", () => hasNoColumn(page, "Age_copy"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A linked column paints its cells with the source column's colours", async () => {
      await session.step(57, "When user colors \"RACE\" column linked to \"AGE\" column", () => colorLinked(page, "RACE", "AGE"));
      await session.step(58, "Then \"RACE\" column should be color-coded linked", () => colorCodedAs(page, "RACE", "linked"));
      await session.step(59, "And the coloring of \"RACE\" column should be linked to \"AGE\" column", () => colorLinkedTo(page, "RACE", "AGE"));
      await session.step(60, "And the \"color of cell 1 of RACE\" reading of grid should be \"#00ff00\"", () => readingReads(page, "color of cell 1 of RACE", el("grid"), "#00ff00"));
      await session.step(61, "And the \"color of cell 3 of RACE\" reading of grid should be \"#ffff00\"", () => readingReads(page, "color of cell 3 of RACE", el("grid"), "#ffff00"));
      await session.step(62, "And the \"color of cell 1 of RACE\" and \"color of cell 3 of RACE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of RACE", "color of cell 3 of RACE", el("grid")));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A change of the source's type leaves the link in place", async () => {
      await session.step(66, "When user colors \"AGE\" column linearly from \"#0000FF\" to \"#FF0000\"", () => colorLinear(page, "AGE", "#0000FF", "#FF0000"));
      await session.step(67, "Then \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
      await session.step(68, "And \"RACE\" column should be color-coded linked", () => colorCodedAs(page, "RACE", "linked"));
      await session.step(69, "And the coloring of \"RACE\" column should be linked to \"AGE\" column", () => colorLinkedTo(page, "RACE", "AGE"));
      await session.step(70, "And the \"color of cell 1 of RACE\" and \"color of cell 1 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of RACE", "color of cell 1 of AGE", el("grid")));
      await session.step(71, "And the \"color of cell 1 of RACE\" reading of grid should not be \"#00ff00\"", () => readingDoesNotRead(page, "color of cell 1 of RACE", el("grid"), "#00ff00"));
      await session.step(72, "When user colors \"AGE\" column conditionally:", () => colorConditional(page, "AGE", [["<40","#00FF00"],[">40","#FF0000"]]));
      await session.step(75, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
      await session.step(76, "And \"RACE\" column should be color-coded linked", () => colorCodedAs(page, "RACE", "linked"));
      await session.step(77, "And the \"color of cell 1 of RACE\" reading of grid should be \"#00ff00\"", () => readingReads(page, "color of cell 1 of RACE", el("grid"), "#00ff00"));
      await session.step(78, "And the \"color of cell 3 of RACE\" reading of grid should be \"#ff0000\"", () => readingReads(page, "color of cell 3 of RACE", el("grid"), "#ff0000"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A linked colouring applied to the text paints the letters, not the cell", async () => {
      await session.step(82, "When user colors the text of \"HEIGHT\" column linked to \"AGE\" column", () => colorLinkedText(page, "HEIGHT", "AGE"));
      await session.step(83, "Then \"HEIGHT\" column should be color-coded linked", () => colorCodedAs(page, "HEIGHT", "linked"));
      await session.step(84, "And the text of \"HEIGHT\" column should be color-coded", () => textColorCoded(page, "HEIGHT"));
      await session.step(85, "And the \"color of cell 1 of HEIGHT\" and \"color of cell 1 of AGE\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of HEIGHT", "color of cell 1 of AGE", el("grid")));
      await session.step(86, "And the \"color of cell 1 of HEIGHT\" and \"color of cell 3 of HEIGHT\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of HEIGHT", "color of cell 3 of HEIGHT", el("grid")));
      await session.step(87, "And the \"cell 1 of HEIGHT\" area of grid should contain the color \"#00FF00\"", () => areaColor(page, "cell 1 of HEIGHT", el("grid"), "#00FF00"));
      await session.step(88, "When user removes the coloring of \"HEIGHT\" column", () => colorOff(page, "HEIGHT"));
      await session.step(89, "And user removes the coloring of \"RACE\" column", () => colorOff(page, "RACE"));
      await session.step(90, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A five-level chain of links all reports Linked", async () => {
      await session.step(94, "When user colors \"SEX\" column linked to \"AGE\" column", () => colorLinked(page, "SEX", "AGE"));
      await session.step(95, "And user colors \"DIS_POP\" column linked to \"SEX\" column", () => colorLinked(page, "DIS_POP", "SEX"));
      await session.step(96, "And user colors \"CONTROL\" column linked to \"DIS_POP\" column", () => colorLinked(page, "CONTROL", "DIS_POP"));
      await session.step(97, "And user colors \"STARTED\" column linked to \"CONTROL\" column", () => colorLinked(page, "STARTED", "CONTROL"));
      await session.step(98, "Then \"SEX\" column should be color-coded linked", () => colorCodedAs(page, "SEX", "linked"));
      await session.step(99, "And \"DIS_POP\" column should be color-coded linked", () => colorCodedAs(page, "DIS_POP", "linked"));
      await session.step(100, "And \"CONTROL\" column should be color-coded linked", () => colorCodedAs(page, "CONTROL", "linked"));
      await session.step(101, "And \"STARTED\" column should be color-coded linked", () => colorCodedAs(page, "STARTED", "linked"));
      await session.step(102, "And the coloring of \"STARTED\" column should be linked to \"CONTROL\" column", () => colorLinkedTo(page, "STARTED", "CONTROL"));
      await session.step(103, "And the \"color of cell 1 of SEX\" reading of grid should be \"#00ff00\"", () => readingReads(page, "color of cell 1 of SEX", el("grid"), "#00ff00"));
      await session.step(104, "And the \"color of cell 3 of SEX\" reading of grid should be \"#ff0000\"", () => readingReads(page, "color of cell 3 of SEX", el("grid"), "#ff0000"));
      await session.step(105, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(106, "And user removes the coloring of \"DIS_POP\" column", () => colorOff(page, "DIS_POP"));
      await session.step(107, "And user removes the coloring of \"CONTROL\" column", () => colorOff(page, "CONTROL"));
      await session.step(108, "And user removes the coloring of \"STARTED\" column", () => colorOff(page, "STARTED"));
      await session.step(109, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Inverting a linear scheme swaps the ends of the gradient", async () => {
      await session.step(112, "When user colors \"AGE\" column linearly from \"#0000FF\" to \"#FF0000\"", () => colorLinear(page, "AGE", "#0000FF", "#FF0000"));
      await session.step(113, "Then the color scheme of \"AGE\" column should be \"#0000FF, #FF0000\"", () => colorSchemeIs(page, "AGE", "#0000FF, #FF0000"));
      await session.step(114, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
      await session.step(115, "When user remembers the \"color of cell 1 of AGE\" reading of grid", () => rememberReading(page, "color of cell 1 of AGE", el("grid")));
      await session.step(116, "And user inverts the color scheme of \"AGE\" column", () => colorInverted(page, "AGE"));
      await session.step(117, "Then the color scheme of \"AGE\" column should be \"#FF0000, #0000FF\"", () => colorSchemeIs(page, "AGE", "#FF0000, #0000FF"));
      await session.step(118, "And grid should have repainted", () => repainted(page, el("grid")));
      await session.step(119, "And the \"color of cell 1 of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "color of cell 1 of AGE", el("grid")));
      await session.step(120, "When user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
      await session.step(121, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(122, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
