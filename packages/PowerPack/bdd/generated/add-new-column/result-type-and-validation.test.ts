/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/result-type-and-validation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.dialogs.add-new-column]
--- */
import {test} from '@playwright/test';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {previewComputes} from '../../bindings/add-new-column.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType, hasColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Add New Column: the result type, a formula the dialog refuses, and a name already taken", () => {
  const session = feature(test, "features/add-new-column/result-type-and-validation.feature", import.meta.url);
  test("Add New Column: the result type, a formula the dialog refuses, and a name already taken", {tag: ["@journey", "@realizes:powerpack.dialogs.add-new-column"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens a table \"calc\" with:", () => openTableOf(page, "calc", [["x","s"],["1","a"],["2","b"],["4","c"]]), [["x","s"],["1","a"],["2","b"],["4","c"]]);
    await session.step(28, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
    await session.step(29, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
    await run.scenario("On auto, the Type choice names the type the formula computes", async () => {
      await session.step(32, "When user types \"${x} / 2\" into formula editor", () => typeInto(page, "${x} / 2", el("formula editor")));
      await session.step(33, "Then the preview grid should show \"${x} / 2\" computed as numbers", () => previewComputes(page, "${x} / 2"));
      await session.step(34, "And column type input should have value \"auto (double)\"", () => shouldHaveValue(page, el("column type input"), "auto (double)"));
      await session.step(35, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(36, "And user types \"${x} > 1\" into formula editor", () => typeInto(page, "${x} > 1", el("formula editor")));
      await session.step(37, "Then column type input should have value \"auto (bool)\"", () => shouldHaveValue(page, el("column type input"), "auto (bool)"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A type picked by hand is the type of the added column", async () => {
      await session.step(41, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(42, "And user types \"dbl\" into column name input", () => typeInto(page, "dbl", el("column name input")));
      await session.step(43, "And user types \"${x} * 2\" into formula editor", () => typeInto(page, "${x} * 2", el("formula editor")));
      await session.step(44, "And user selects \"double\" in column type input", () => selectIn(page, "double", el("column type input")));
      await session.step(45, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(46, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(47, "And \"dbl\" column should have type \"double\"", () => columnType(page, "dbl", "double"));
      await session.step(48, "And the value of \"dbl\" column in row 3 should be \"8\"", () => valueInRow(page, "dbl", 3, "8"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A picked type the formula cannot give is named, and OK is disabled", async () => {
      await session.step(52, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(53, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(54, "When user types \"${s}\" into formula editor", () => typeInto(page, "${s}", el("formula editor")));
      await session.step(55, "And user selects \"int\" in column type input", () => selectIn(page, "int", el("column type input")));
      await session.step(56, "Then formula error should contain text \"doesn't match with current column type (int)\"", () => shouldContainText(page, el("formula error"), "doesn't match with current column type (int)"));
      await session.step(57, "And OK button in \"Add New Column\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add New Column\" dialog"), "disabled"));
      await session.step(58, "When user selects \"string\" in column type input", () => selectIn(page, "string", el("column type input")));
      await session.step(59, "Then OK button in \"Add New Column\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Add New Column\" dialog"), "enabled"));
      await session.step(60, "When user clicks on CANCEL button in \"Add New Column\" dialog", () => clickOn(page, el("CANCEL button in \"Add New Column\" dialog")));
      await session.step(61, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A formula naming a column the table has not got is refused", async () => {
      await session.step(65, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(66, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(67, "When user types \"${nope} + 1\" into formula editor", () => typeInto(page, "${nope} + 1", el("formula editor")));
      await session.step(68, "Then formula error should contain text \"nope\"", () => shouldContainText(page, el("formula error"), "nope"));
      await session.step(69, "And OK button in \"Add New Column\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Add New Column\" dialog"), "disabled"));
      await session.step(70, "When user clicks on CANCEL button in \"Add New Column\" dialog", () => clickOn(page, el("CANCEL button in \"Add New Column\" dialog")));
      await session.step(71, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A vector function inside arithmetic is a valid formula and adds its column", async () => {
      await session.step(75, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(76, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(77, "When user types \"running\" into column name input", () => typeInto(page, "running", el("column name input")));
      await session.step(78, "And user types \"CumSum(${x}) * 2\" into formula editor", () => typeInto(page, "CumSum(${x}) * 2", el("formula editor")));
      await session.step(79, "Then the preview grid should show \"running\" computed as numbers", () => previewComputes(page, "running"));
      await session.step(80, "When user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(81, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(82, "And the value of \"running\" column in row 1 should be \"2\"", () => valueInRow(page, "running", 1, "2"));
      await session.step(83, "And the value of \"running\" column in row 3 should be \"14\"", () => valueInRow(page, "running", 3, "14"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Plain text makes the formula a text template filled in for every row", async () => {
      await session.step(87, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(88, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(89, "When user types \"template\" into column name input", () => typeInto(page, "template", el("column name input")));
      await session.step(90, "And user types \"${s} + ${s}\" into formula editor", () => typeInto(page, "${s} + ${s}", el("formula editor")));
      await session.step(91, "And user selects \"plain text\" in column type input", () => selectIn(page, "plain text", el("column type input")));
      await session.step(92, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(93, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(94, "And \"template\" column should have type \"string\"", () => columnType(page, "template", "string"));
      await session.step(95, "And the value of \"template\" column in row 1 should be \"a + a\"", () => valueInRow(page, "template", 1, "a + a"));
      await session.step(96, "And the value of \"template\" column in row 3 should be \"c + c\"", () => valueInRow(page, "template", 3, "c + c"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A name the table already has is kept apart, and OK adds it under the new name", async () => {
      await session.step(100, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(101, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(102, "When user types \"${x} + 10\" into formula editor", () => typeInto(page, "${x} + 10", el("formula editor")));
      await session.step(103, "Then the preview grid should show \"${x} + 10\" computed as numbers", () => previewComputes(page, "${x} + 10"));
      await session.step(104, "When user types \"x\" into column name input", () => typeInto(page, "x", el("column name input")));
      await session.step(105, "Then formula error should contain text \"Column named x already exists. Column will be named x (2)\"", () => shouldContainText(page, el("formula error"), "Column named x already exists. Column will be named x (2)"));
      await session.step(106, "When user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(107, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(108, "And the table should have a column \"x (2)\"", () => hasColumn(page, "x (2)"));
      await session.step(109, "And the value of \"x (2)\" column in row 3 should be \"14\"", () => valueInRow(page, "x (2)", 3, "14"));
      await session.step(110, "And the value of \"x\" column in row 3 should be \"4\"", () => valueInRow(page, "x", 3, "4"));
      await session.step(111, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A number the formula computes is added as text with string, and as the formula's text with plain text", async () => {
      await session.step(114, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(115, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(116, "When user types \"as text\" into column name input", () => typeInto(page, "as text", el("column name input")));
      await session.step(117, "And user selects \"string\" in column type input", () => selectIn(page, "string", el("column type input")));
      await session.step(118, "And user types \"${x} + 1\" into formula editor", () => typeInto(page, "${x} + 1", el("formula editor")));
      await session.step(119, "Then the preview grid should show \"as text\" computed as numbers", () => previewComputes(page, "as text"));
      await session.step(120, "When user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(121, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(122, "And \"as text\" column should have type \"string\"", () => columnType(page, "as text", "string"));
      await session.step(123, "And the value of \"as text\" column in row 3 should be \"5\"", () => valueInRow(page, "as text", 3, "5"));
      await session.step(124, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(125, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(126, "When user types \"as template\" into column name input", () => typeInto(page, "as template", el("column name input")));
      await session.step(127, "And user selects \"plain text\" in column type input", () => selectIn(page, "plain text", el("column type input")));
      await session.step(128, "And user types \"${x} + 1\" into formula editor", () => typeInto(page, "${x} + 1", el("formula editor")));
      await session.step(129, "Then the \"text of cell 3 of as template\" reading of preview grid viewer should be \"4 + 1\"", () => readingReads(page, "text of cell 3 of as template", el("preview grid viewer"), "4 + 1"));
      await session.step(130, "When user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(131, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(132, "And \"as template\" column should have type \"string\"", () => columnType(page, "as template", "string"));
      await session.step(133, "And the value of \"as template\" column in row 3 should be \"4 + 1\"", () => valueInRow(page, "as template", 3, "4 + 1"));
      await session.step(134, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With the formula left empty, OK adds an empty column of the picked type", async () => {
      await session.step(137, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(138, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(139, "When user selects \"int\" in column type input", () => selectIn(page, "int", el("column type input")));
      await session.step(140, "And user types \"notes\" into column name input", () => typeInto(page, "notes", el("column name input")));
      await session.step(141, "Then the \"text of cell 1 of notes\" reading of preview grid viewer should be \"\"", () => readingReads(page, "text of cell 1 of notes", el("preview grid viewer"), ""));
      await session.step(142, "When user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(143, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(144, "And \"notes\" column should have type \"int\"", () => columnType(page, "notes", "int"));
      await session.step(145, "And the value of \"notes\" column in row 1 should be \"\"", () => valueInRow(page, "notes", 1, ""));
      await session.step(146, "And the value of \"notes\" column in row 3 should be \"\"", () => valueInRow(page, "notes", 3, ""));
      await session.step(147, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
