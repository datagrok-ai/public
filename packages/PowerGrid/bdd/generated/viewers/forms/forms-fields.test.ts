/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-fields.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {renameColumn} from '../../../bindings/forms.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn, hasNoColumn, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hasArea, hasNoArea, noBalloons, noErrors, propertyShouldBe, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer field set and number format", () => {
  const session = feature(test, "features/viewers/forms/forms-fields.feature", import.meta.url);
  test("Forms viewer field set and number format", {tag: ["@journey", "@viewers", "@realizes:viewers.forms"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a forms viewer with:", () => addViewerWith(page, "forms", [["Show Mouse Over Row","false"]]));
    await session.step(14, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
    await session.step(15, "And the \"fields shown\" reading of forms viewer should be 11", () => readingIs(page, "fields shown", el("forms viewer"), 11));
    await run.scenario("The fields are the picked columns, in the picked order", async () => {
      await session.step(18, "When user sets \"Fields\" property of forms viewer to \"RACE, AGE, SEX\"", () => setProperty(page, "Fields", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(19, "Then the \"fields\" reading of forms viewer should be \"RACE, AGE, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(20, "And the \"fields shown\" reading of forms viewer should be 3", () => readingIs(page, "fields shown", el("forms viewer"), 3));
      await session.step(21, "And \"Fields\" property of forms viewer should be \"RACE, AGE, SEX\"", () => propertyShouldBe(page, "Fields", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(22, "And forms viewer should have a \"label RACE\" area", () => hasArea(page, el("forms viewer"), "label RACE"));
      await session.step(23, "And forms viewer should have a \"label AGE\" area", () => hasArea(page, el("forms viewer"), "label AGE"));
      await session.step(24, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The header's remove icon drops a field", async () => {
      await session.step(27, "When user clicks on the \"remove AGE\" area of forms viewer", () => clickArea(page, "remove AGE", el("forms viewer")));
      await session.step(28, "Then the \"fields\" reading of forms viewer should be \"RACE, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "RACE, SEX"));
      await session.step(29, "And the \"fields shown\" reading of forms viewer should be 2", () => readingIs(page, "fields shown", el("forms viewer"), 2));
      await session.step(30, "And forms viewer should not have a \"label AGE\" area", () => hasNoArea(page, el("forms viewer"), "label AGE"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column renamed to a \"~\" name leaves the field set", async () => {
      await session.step(34, "When user renames \"RACE\" column to \"~RACE\"", () => renameColumn(page, "RACE", "~RACE"));
      await session.step(35, "Then the \"fields\" reading of forms viewer should be \"SEX\"", () => readingReads(page, "fields", el("forms viewer"), "SEX"));
      await session.step(36, "And the \"fields shown\" reading of forms viewer should be 1", () => readingIs(page, "fields shown", el("forms viewer"), 1));
      await session.step(37, "And forms viewer should not have a \"label RACE\" area", () => hasNoArea(page, el("forms viewer"), "label RACE"));
      await session.step(38, "When user renames \"~RACE\" column to \"RACE\"", () => renameColumn(page, "~RACE", "RACE"));
      await session.step(39, "Then the table should have a column \"RACE\"", () => hasColumn(page, "RACE"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A named number format reaches the float fields only", async () => {
      await session.step(43, "When user sets \"Fields\" property of forms viewer to \"HEIGHT, AGE, SEX\"", () => setProperty(page, "Fields", el("forms viewer"), "HEIGHT, AGE, SEX"));
      await session.step(44, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(45, "Then \"Number Format\" property of forms viewer should be \"Same as grid\"", () => propertyShouldBe(page, "Number Format", el("forms viewer"), "Same as grid"));
      await session.step(46, "And the \"AGE of current card\" reading of forms viewer should be \"26\"", () => readingReads(page, "AGE of current card", el("forms viewer"), "26"));
      await session.step(47, "And the \"SEX of current card\" reading of forms viewer should be \"F\"", () => readingReads(page, "SEX of current card", el("forms viewer"), "F"));
      await session.step(48, "When user sets \"Number Format\" property of forms viewer to \"3 significant digits\"", () => setProperty(page, "Number Format", el("forms viewer"), "3 significant digits"));
      await session.step(49, "Then the \"HEIGHT of current card\" reading of forms viewer should be \"175\"", () => readingReads(page, "HEIGHT of current card", el("forms viewer"), "175"));
      await session.step(50, "And the \"AGE of current card\" reading of forms viewer should be \"26\"", () => readingReads(page, "AGE of current card", el("forms viewer"), "26"));
      await session.step(51, "And the \"SEX of current card\" reading of forms viewer should be \"F\"", () => readingReads(page, "SEX of current card", el("forms viewer"), "F"));
      await session.step(52, "When user sets \"Number Format\" property of forms viewer to \"3 digits after comma\"", () => setProperty(page, "Number Format", el("forms viewer"), "3 digits after comma"));
      await session.step(53, "Then the \"HEIGHT of current card\" reading of forms viewer should be \"174.705\"", () => readingReads(page, "HEIGHT of current card", el("forms viewer"), "174.705"));
      await session.step(54, "When user sets \"Number Format\" property of forms viewer to \"Same as grid\"", () => setProperty(page, "Number Format", el("forms viewer"), "Same as grid"));
      await session.step(55, "Then \"Number Format\" property of forms viewer should be \"Same as grid\"", () => propertyShouldBe(page, "Number Format", el("forms viewer"), "Same as grid"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An empty field set draws no card and says nothing", async () => {
      await session.step(59, "When user sets \"Fields\" property of forms viewer to \"\"", () => setProperty(page, "Fields", el("forms viewer"), ""));
      await session.step(60, "Then the \"fields shown\" reading of forms viewer should be 0", () => readingIs(page, "fields shown", el("forms viewer"), 0));
      await session.step(61, "And the \"cards\" reading of forms viewer should be 0", () => readingIs(page, "cards", el("forms viewer"), 0));
      await session.step(62, "And forms viewer should not have a \"label HEIGHT\" area", () => hasNoArea(page, el("forms viewer"), "label HEIGHT"));
      await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(64, "When user sets \"Fields\" property of forms viewer to \"USUBJID, AGE, SEX\"", () => setProperty(page, "Fields", el("forms viewer"), "USUBJID, AGE, SEX"));
      await session.step(65, "Then the \"fields shown\" reading of forms viewer should be 3", () => readingIs(page, "fields shown", el("forms viewer"), 3));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column removed from the table takes its field with it", async () => {
      await session.step(69, "When user removes \"AGE\" column", () => removeColumn(page, "AGE"));
      await session.step(70, "Then the table should not have a column \"AGE\"", () => hasNoColumn(page, "AGE"));
      await session.step(71, "And the \"fields\" reading of forms viewer should be \"USUBJID, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "USUBJID, SEX"));
      await session.step(72, "And forms viewer should not have a \"label AGE\" area", () => hasNoArea(page, el("forms viewer"), "label AGE"));
      await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
