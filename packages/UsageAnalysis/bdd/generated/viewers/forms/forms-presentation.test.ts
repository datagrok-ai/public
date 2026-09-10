/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-presentation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorCodedAs, colorLinear, colorOff, noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, propertyShouldBe, readingAsRemembered, readingDoesNotRead, readingNotAsRemembered, readingReads, readingsDiffer, readingsEqual, rememberReading, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer colour coding, alignment and font", () => {
  const session = feature(test, "features/viewers/forms/forms-presentation.feature", import.meta.url);
  test("Forms viewer colour coding, alignment and font", {tag: ["@journey", "@viewers", "@realizes:viewers.forms"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user colors \"AGE\" column linearly from \"#FF0000\" to \"#00FF00\"", () => colorLinear(page, "AGE", "#FF0000", "#00FF00"));
    await session.step(19, "And user adds a forms viewer", () => addViewer(page, "forms"));
    await session.step(20, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(21, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
    await session.step(22, "And \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
    await session.step(23, "And the \"AGE of card 1\" reading of forms viewer should be \"26\"", () => readingReads(page, "AGE of card 1", el("forms viewer"), "26"));
    await run.scenario("A colour-coded column paints the field's background, an uncoded one does not", async () => {
      await session.step(26, "Then \"colorCode\" property of forms viewer should be \"true\"", () => propertyShouldBe(page, "colorCode", el("forms viewer"), "true"));
      await session.step(27, "And the \"background of USUBJID of card 1\" reading of forms viewer should be \"#FFFFFF\"", () => readingReads(page, "background of USUBJID of card 1", el("forms viewer"), "#FFFFFF"));
      await session.step(28, "And the \"background of AGE of card 1\" reading of forms viewer should not be \"#FFFFFF\"", () => readingDoesNotRead(page, "background of AGE of card 1", el("forms viewer"), "#FFFFFF"));
      await session.step(29, "And the \"background of AGE of card 1\" and \"background of USUBJID of card 1\" readings of forms viewer should differ", () => readingsDiffer(page, "background of AGE of card 1", "background of USUBJID of card 1", el("forms viewer")));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color Code off drops the background to the uncoded one, and on restores it", async () => {
      await session.step(33, "When user sets \"colorCode\" property of forms viewer to \"false\"", () => setProperty(page, "colorCode", el("forms viewer"), "false"));
      await session.step(34, "Then the \"background of AGE of card 1\" and \"background of USUBJID of card 1\" readings of forms viewer should be the same", () => readingsEqual(page, "background of AGE of card 1", "background of USUBJID of card 1", el("forms viewer")));
      await session.step(35, "And the \"background of AGE of card 1\" reading of forms viewer should be \"#FFFFFF\"", () => readingReads(page, "background of AGE of card 1", el("forms viewer"), "#FFFFFF"));
      await session.step(36, "When user sets \"colorCode\" property of forms viewer to \"true\"", () => setProperty(page, "colorCode", el("forms viewer"), "true"));
      await session.step(37, "Then the \"background of AGE of card 1\" and \"background of USUBJID of card 1\" readings of forms viewer should differ", () => readingsDiffer(page, "background of AGE of card 1", "background of USUBJID of card 1", el("forms viewer")));
      await session.step(38, "And the \"background of AGE of card 1\" reading of forms viewer should not be \"#FFFFFF\"", () => readingDoesNotRead(page, "background of AGE of card 1", el("forms viewer"), "#FFFFFF"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A row of another colour gives its field another background", async () => {
      await session.step(42, "Given user remembers the \"background of AGE of card 1\" reading of forms viewer", () => rememberReading(page, "background of AGE of card 1", el("forms viewer")));
      await session.step(43, "When user makes row 78 current", () => makeRowCurrent(page, 78));
      await session.step(44, "Then the \"AGE of card 1\" reading of forms viewer should be \"60\"", () => readingReads(page, "AGE of card 1", el("forms viewer"), "60"));
      await session.step(45, "And the \"background of AGE of card 1\" reading of forms viewer should not be as remembered", () => readingNotAsRemembered(page, "background of AGE of card 1", el("forms viewer")));
      await session.step(46, "And the \"background of AGE of card 1\" reading of forms viewer should not be \"#FFFFFF\"", () => readingDoesNotRead(page, "background of AGE of card 1", el("forms viewer"), "#FFFFFF"));
      await session.step(47, "When user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(48, "Then the \"background of AGE of card 1\" reading of forms viewer should be as remembered", () => readingAsRemembered(page, "background of AGE of card 1", el("forms viewer")));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the colouring leaves every field the same background", async () => {
      await session.step(52, "When user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
      await session.step(53, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(54, "And the \"background of AGE of card 1\" and \"background of USUBJID of card 1\" readings of forms viewer should be the same", () => readingsEqual(page, "background of AGE of card 1", "background of USUBJID of card 1", el("forms viewer")));
      await session.step(55, "When user colors \"AGE\" column linearly from \"#FF0000\" to \"#00FF00\"", () => colorLinear(page, "AGE", "#FF0000", "#00FF00"));
      await session.step(56, "Then the \"background of AGE of card 1\" and \"background of USUBJID of card 1\" readings of forms viewer should differ", () => readingsDiffer(page, "background of AGE of card 1", "background of USUBJID of card 1", el("forms viewer")));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Font property reaches every field, and the alignment is the stylesheet's", async () => {
      await session.step(60, "Then the \"align of AGE of card 1\" reading of forms viewer should be \"start\"", () => readingReads(page, "align of AGE of card 1", el("forms viewer"), "start"));
      await session.step(61, "And the \"align of SEX of card 1\" reading of forms viewer should be \"start\"", () => readingReads(page, "align of SEX of card 1", el("forms viewer"), "start"));
      await session.step(62, "And the \"font of AGE of card 1\" reading of forms viewer should be \"13px Roboto\"", () => readingReads(page, "font of AGE of card 1", el("forms viewer"), "13px Roboto"));
      await session.step(63, "When user sets \"font\" property of forms viewer to \"italic bold 15px \\\"Times New Roman\\\"\"", () => setProperty(page, "font", el("forms viewer"), "italic bold 15px \"Times New Roman\""));
      await session.step(64, "Then the \"font of AGE of card 1\" reading of forms viewer should be \"italic bold 15px \\\"Times New Roman\\\"\"", () => readingReads(page, "font of AGE of card 1", el("forms viewer"), "italic bold 15px \"Times New Roman\""));
      await session.step(65, "And the \"font of SEX of card 1\" reading of forms viewer should be \"italic bold 15px \\\"Times New Roman\\\"\"", () => readingReads(page, "font of SEX of card 1", el("forms viewer"), "italic bold 15px \"Times New Roman\""));
      await session.step(66, "And the \"align of AGE of card 1\" reading of forms viewer should be \"start\"", () => readingReads(page, "align of AGE of card 1", el("forms viewer"), "start"));
      await session.step(67, "When user sets \"font\" property of forms viewer to \"normal normal 13px \\\"Roboto\\\"\"", () => setProperty(page, "font", el("forms viewer"), "normal normal 13px \"Roboto\""));
      await session.step(68, "Then the \"font of AGE of card 1\" reading of forms viewer should be \"13px Roboto\"", () => readingReads(page, "font of AGE of card 1", el("forms viewer"), "13px Roboto"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
