/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-fields.feature
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
import {columnCount, displayedInRow, hasColumn, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, removeColumn, renameColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, hasArea, hasNoArea, noBalloons, noErrors, propertyShouldBe, readingIs, readingReads, readingsEqual, reportsNoError, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer field lifecycle and number format", () => {
  const session = feature(test, "features/viewers/forms/forms-fields.feature", import.meta.url);
  test("Forms viewer field lifecycle and number format", {tag: ["@journey", "@viewers", "@realizes:viewers.forms"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(27, "And user adds a calculated column \"COMPUTED_H\" with formula \"${HEIGHT}\"", () => addCalculated(page, "COMPUTED_H", "${HEIGHT}"));
    await session.step(28, "And user adds a forms viewer", () => addViewer(page, "forms"));
    await session.step(29, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
    await session.step(30, "And the table should have 12 columns", () => columnCount(page, 12));
    await session.step(31, "And the \"fields shown\" reading of forms viewer should be 12", () => readingIs(page, "fields shown", el("forms viewer"), 12));
    await session.step(32, "And the \"fields\" and \"header labels\" readings of forms viewer should be the same", () => readingsEqual(page, "fields", "header labels", el("forms viewer")));
    await run.scenario("The fields are drawn in the order they were given, not in table order", async () => {
      await session.step(35, "When user sets \"fieldsColumnNames\" property of forms viewer to \"RACE, AGE, SEX\"", () => setProperty(page, "fieldsColumnNames", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(36, "Then the \"fields\" reading of forms viewer should be \"RACE, AGE, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(37, "And the \"header labels\" reading of forms viewer should be \"RACE, AGE, SEX\"", () => readingReads(page, "header labels", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(38, "And the \"fields shown\" reading of forms viewer should be 3", () => readingIs(page, "fields shown", el("forms viewer"), 3));
      await session.step(39, "And forms viewer should have a \"field RACE of card 1\" area", () => hasArea(page, el("forms viewer"), "field RACE of card 1"));
      await session.step(40, "And forms viewer should not have a \"field WEIGHT of card 1\" area", () => hasNoArea(page, el("forms viewer"), "field WEIGHT of card 1"));
      await session.step(41, "And the \"RACE of card 1\" reading of forms viewer should be \"Caucasian\"", () => readingReads(page, "RACE of card 1", el("forms viewer"), "Caucasian"));
      await session.step(42, "When user sets \"fieldsColumnNames\" property of forms viewer to \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, COMPUTED_H\"", () => setProperty(page, "fieldsColumnNames", el("forms viewer"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, COMPUTED_H"));
      await session.step(43, "Then the \"fields shown\" reading of forms viewer should be 12", () => readingIs(page, "fields shown", el("forms viewer"), 12));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The header cross takes a field out and leaves the order of the rest", async () => {
      await session.step(47, "When user sets \"fieldsColumnNames\" property of forms viewer to \"RACE, AGE, SEX\"", () => setProperty(page, "fieldsColumnNames", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(48, "Then forms viewer should have a \"remove AGE\" area", () => hasArea(page, el("forms viewer"), "remove AGE"));
      await session.step(49, "When user clicks on the \"remove AGE\" area of forms viewer", () => clickArea(page, "remove AGE", el("forms viewer")));
      await session.step(50, "Then the \"fields\" reading of forms viewer should be \"RACE, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "RACE, SEX"));
      await session.step(51, "And the \"header labels\" reading of forms viewer should be \"RACE, SEX\"", () => readingReads(page, "header labels", el("forms viewer"), "RACE, SEX"));
      await session.step(52, "And the \"fields shown\" reading of forms viewer should be 2", () => readingIs(page, "fields shown", el("forms viewer"), 2));
      await session.step(53, "And forms viewer should not have a \"field AGE of card 1\" area", () => hasNoArea(page, el("forms viewer"), "field AGE of card 1"));
      await session.step(54, "And forms viewer should not have a \"label AGE\" area", () => hasNoArea(page, el("forms viewer"), "label AGE"));
      await session.step(55, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Number Format is the viewer's own and leaves the value alone", async () => {
      await session.step(59, "When user sets \"fieldsColumnNames\" property of forms viewer to \"COMPUTED_H, AGE, SEX\"", () => setProperty(page, "fieldsColumnNames", el("forms viewer"), "COMPUTED_H, AGE, SEX"));
      await session.step(60, "Then \"numberFormat\" property of forms viewer should be \"Same as grid\"", () => propertyShouldBe(page, "numberFormat", el("forms viewer"), "Same as grid"));
      await session.step(61, "And the \"COMPUTED_H of card 1\" reading of forms viewer should be \"174.71\"", () => readingReads(page, "COMPUTED_H of card 1", el("forms viewer"), "174.71"));
      await session.step(62, "And the \"COMPUTED_H\" cell of row 1 should be displayed as \"174.71\"", () => displayedInRow(page, "COMPUTED_H", 1, "174.71"));
      await session.step(63, "And the \"AGE of card 1\" reading of forms viewer should be \"26\"", () => readingReads(page, "AGE of card 1", el("forms viewer"), "26"));
      await session.step(64, "And the \"SEX of card 1\" reading of forms viewer should be \"F\"", () => readingReads(page, "SEX of card 1", el("forms viewer"), "F"));
      await session.step(65, "When user sets \"numberFormat\" property of forms viewer to \"3 digits after comma\"", () => setProperty(page, "numberFormat", el("forms viewer"), "3 digits after comma"));
      await session.step(66, "Then the \"COMPUTED_H of card 1\" reading of forms viewer should be \"174.705\"", () => readingReads(page, "COMPUTED_H of card 1", el("forms viewer"), "174.705"));
      await session.step(67, "And the \"COMPUTED_H\" cell of row 1 should be displayed as \"174.71\"", () => displayedInRow(page, "COMPUTED_H", 1, "174.71"));
      await session.step(68, "And the \"AGE of card 1\" reading of forms viewer should be \"26\"", () => readingReads(page, "AGE of card 1", el("forms viewer"), "26"));
      await session.step(69, "And the \"SEX of card 1\" reading of forms viewer should be \"F\"", () => readingReads(page, "SEX of card 1", el("forms viewer"), "F"));
      await session.step(70, "And the value of \"COMPUTED_H\" column in row 1 should be \"174.7050018310547\"", () => valueInRow(page, "COMPUTED_H", 1, "174.7050018310547"));
      await session.step(71, "When user sets \"numberFormat\" property of forms viewer to \"2 digits after comma\"", () => setProperty(page, "numberFormat", el("forms viewer"), "2 digits after comma"));
      await session.step(72, "Then the \"COMPUTED_H of card 1\" reading of forms viewer should be \"174.71\"", () => readingReads(page, "COMPUTED_H of card 1", el("forms viewer"), "174.71"));
      await session.step(73, "When user sets \"numberFormat\" property of forms viewer to \"Same as grid\"", () => setProperty(page, "numberFormat", el("forms viewer"), "Same as grid"));
      await session.step(74, "Then the \"COMPUTED_H of card 1\" reading of forms viewer should be \"174.71\"", () => readingReads(page, "COMPUTED_H of card 1", el("forms viewer"), "174.71"));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An empty field set draws nothing and says nothing", async () => {
      await session.step(78, "When user sets \"fieldsColumnNames\" property of forms viewer to \"\"", () => setProperty(page, "fieldsColumnNames", el("forms viewer"), ""));
      await session.step(79, "Then the \"fields shown\" reading of forms viewer should be 0", () => readingIs(page, "fields shown", el("forms viewer"), 0));
      await session.step(80, "And the \"fields\" reading of forms viewer should be \"\"", () => readingReads(page, "fields", el("forms viewer"), ""));
      await session.step(81, "And the \"header labels\" reading of forms viewer should be \"\"", () => readingReads(page, "header labels", el("forms viewer"), ""));
      await session.step(82, "And forms viewer should not have a \"label AGE\" area", () => hasNoArea(page, el("forms viewer"), "label AGE"));
      await session.step(83, "And forms viewer should report no error", () => reportsNoError(page, el("forms viewer")));
      await session.step(84, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
      await session.step(86, "When user sets \"fieldsColumnNames\" property of forms viewer to \"RACE, AGE, SEX\"", () => setProperty(page, "fieldsColumnNames", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(87, "Then the \"fields shown\" reading of forms viewer should be 3", () => readingIs(page, "fields shown", el("forms viewer"), 3));
    });
    await run.scenario("A renamed column carries its field to the new name", async () => {
      await session.step(90, "Then the \"fields\" reading of forms viewer should be \"RACE, AGE, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(91, "When user renames \"SEX\" column to \"GENDER\"", () => renameColumn(page, "SEX", "GENDER"));
      await session.step(92, "Then the table should have a column \"GENDER\"", () => hasColumn(page, "GENDER"));
      await session.step(93, "And the \"fields\" reading of forms viewer should be \"RACE, AGE, GENDER\"", () => readingReads(page, "fields", el("forms viewer"), "RACE, AGE, GENDER"));
      await session.step(94, "And the \"header labels\" reading of forms viewer should be \"RACE, AGE, GENDER\"", () => readingReads(page, "header labels", el("forms viewer"), "RACE, AGE, GENDER"));
      await session.step(95, "And forms viewer should have a \"field GENDER of card 1\" area", () => hasArea(page, el("forms viewer"), "field GENDER of card 1"));
      await session.step(96, "And forms viewer should not have a \"field SEX of card 1\" area", () => hasNoArea(page, el("forms viewer"), "field SEX of card 1"));
      await session.step(97, "And the \"GENDER of card 1\" reading of forms viewer should be \"F\"", () => readingReads(page, "GENDER of card 1", el("forms viewer"), "F"));
      await session.step(98, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(99, "When user renames \"GENDER\" column to \"SEX\"", () => renameColumn(page, "GENDER", "SEX"));
      await session.step(100, "Then the \"fields\" reading of forms viewer should be \"RACE, AGE, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A rename to a service name drops the field", async () => {
      await session.step(104, "Then the \"fields\" reading of forms viewer should be \"RACE, AGE, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(105, "When user renames \"SEX\" column to \"~SERVICE\"", () => renameColumn(page, "SEX", "~SERVICE"));
      await session.step(106, "Then the \"fields\" reading of forms viewer should be \"RACE, AGE\"", () => readingReads(page, "fields", el("forms viewer"), "RACE, AGE"));
      await session.step(107, "And the \"header labels\" reading of forms viewer should be \"RACE, AGE\"", () => readingReads(page, "header labels", el("forms viewer"), "RACE, AGE"));
      await session.step(108, "And the \"header labels\" reading of forms viewer should not contain \"~\"", () => readingNotContains(page, "header labels", el("forms viewer"), "~"));
      await session.step(109, "And the \"fields shown\" reading of forms viewer should be 2", () => readingIs(page, "fields shown", el("forms viewer"), 2));
      await session.step(110, "And forms viewer should not have a \"field ~SERVICE of card 1\" area", () => hasNoArea(page, el("forms viewer"), "field ~SERVICE of card 1"));
      await session.step(111, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(112, "When user renames \"~SERVICE\" column to \"SEX\"", () => renameColumn(page, "~SERVICE", "SEX"));
      await session.step(113, "Then the table should have a column \"SEX\"", () => hasColumn(page, "SEX"));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column that leaves the table takes its field with it", async () => {
      await session.step(117, "When user sets \"fieldsColumnNames\" property of forms viewer to \"RACE, AGE, SEX\"", () => setProperty(page, "fieldsColumnNames", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(118, "Then the \"fields\" reading of forms viewer should be \"RACE, AGE, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "RACE, AGE, SEX"));
      await session.step(119, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
      await session.step(120, "Then the table should not have a column \"RACE\"", () => hasNoColumn(page, "RACE"));
      await session.step(121, "And the \"fields\" reading of forms viewer should be \"AGE, SEX\"", () => readingReads(page, "fields", el("forms viewer"), "AGE, SEX"));
      await session.step(122, "And the \"header labels\" reading of forms viewer should be \"AGE, SEX\"", () => readingReads(page, "header labels", el("forms viewer"), "AGE, SEX"));
      await session.step(123, "And the \"fields shown\" reading of forms viewer should be 2", () => readingIs(page, "fields shown", el("forms viewer"), 2));
      await session.step(124, "And forms viewer should not have a \"field RACE of card 1\" area", () => hasNoArea(page, el("forms viewer"), "field RACE of card 1"));
      await session.step(125, "And the \"AGE of card 1\" reading of forms viewer should be \"26\"", () => readingReads(page, "AGE of card 1", el("forms viewer"), "26"));
      await session.step(126, "And forms viewer should report no error", () => reportsNoError(page, el("forms viewer")));
      await session.step(127, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(128, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
