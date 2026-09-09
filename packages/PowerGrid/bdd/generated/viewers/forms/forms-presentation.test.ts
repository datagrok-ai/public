/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-presentation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clearSelection, colorLinear, colorOff, selectWhereOneOf} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, readingAsRemembered, readingDiffers, readingIs, readingReads, rememberReading, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer presentation", () => {
  const session = feature(test, "features/viewers/forms/forms-presentation.feature", import.meta.url);
  test("Forms viewer presentation", {tag: ["@journey", "@viewers", "@realizes:viewers.forms"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a forms viewer with:", () => addViewerWith(page, "forms", [["Show Mouse Over Row","false"]]));
    await session.step(14, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(15, "Then the \"USUBJID of current card\" reading of forms viewer should be \"X0273T21000300003\"", () => readingReads(page, "USUBJID of current card", el("forms viewer"), "X0273T21000300003"));
    await run.scenario("A text column is an input field", async () => {
      await session.step(18, "Then the \"field kind of USUBJID\" reading of forms viewer should be \"input\"", () => readingReads(page, "field kind of USUBJID", el("forms viewer"), "input"));
      await session.step(19, "And the \"field kind of AGE\" reading of forms viewer should be \"input\"", () => readingReads(page, "field kind of AGE", el("forms viewer"), "input"));
      await session.step(20, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column's colour coding paints its field and Color Code gates it", async () => {
      await session.step(23, "When user colors \"AGE\" column linearly from \"#FF0000\" to \"#0000FF\"", () => colorLinear(page, "AGE", "#FF0000", "#0000FF"));
      await session.step(24, "Then the \"background of AGE of current card\" reading of forms viewer should differ from before", () => readingDiffers(page, "background of AGE of current card", el("forms viewer")));
      await session.step(25, "When user remembers the \"background of AGE of current card\" reading of forms viewer", () => rememberReading(page, "background of AGE of current card", el("forms viewer")));
      await session.step(26, "And user sets \"Color Code\" property of forms viewer to \"false\"", () => setProperty(page, "Color Code", el("forms viewer"), "false"));
      await session.step(27, "Then the \"background of AGE of current card\" reading of forms viewer should differ from before", () => readingDiffers(page, "background of AGE of current card", el("forms viewer")));
      await session.step(28, "When user sets \"Color Code\" property of forms viewer to \"true\"", () => setProperty(page, "Color Code", el("forms viewer"), "true"));
      await session.step(29, "Then the \"background of AGE of current card\" reading of forms viewer should be as remembered", () => readingAsRemembered(page, "background of AGE of current card", el("forms viewer")));
      await session.step(30, "When user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
      await session.step(31, "Then the \"background of AGE of current card\" reading of forms viewer should differ from before", () => readingDiffers(page, "background of AGE of current card", el("forms viewer")));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The colouring reaches the selected rows' cards too", async () => {
      await session.step(35, "When user selects rows where \"USUBJID\" is one of \"X0273T21000400001, X0273T21001500015\"", () => selectWhereOneOf(page, "USUBJID", "X0273T21000400001, X0273T21001500015"));
      await session.step(36, "Then the \"cards\" reading of forms viewer should be 3", () => readingIs(page, "cards", el("forms viewer"), 3));
      await session.step(37, "When user colors \"AGE\" column linearly from \"#FF0000\" to \"#0000FF\"", () => colorLinear(page, "AGE", "#FF0000", "#0000FF"));
      await session.step(38, "Then the \"background of AGE of card 2\" reading of forms viewer should differ from before", () => readingDiffers(page, "background of AGE of card 2", el("forms viewer")));
      await session.step(39, "And the \"background of AGE of card 3\" reading of forms viewer should differ from before", () => readingDiffers(page, "background of AGE of card 3", el("forms viewer")));
      await session.step(40, "When user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
      await session.step(41, "And user clears the row selection", () => clearSelection(page));
      await session.step(42, "Then the \"cards\" reading of forms viewer should be 1", () => readingIs(page, "cards", el("forms viewer"), 1));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
