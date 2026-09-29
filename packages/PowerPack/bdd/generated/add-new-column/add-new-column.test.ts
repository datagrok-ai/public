/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/add-new-column.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.dialogs.add-new-column, powerpack.cp.add-new-column-persists]
--- */
import {test} from '@playwright/test';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '../../bindings/io.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {acceptCompletion, dragAreaOnto, dragColumnName, dragCornerBy, holdsFormula, keepsWidth, previewComputes, rememberSize, sizeAgainstRemembered, typeAtCaret} from '../../bindings/add-new-column.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, pressKeyIn, shouldBe, shouldContainText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {liesWithin} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Add New Column on demog: the dialog, a formula built by hand, and the input history", () => {
  const session = feature(test, "features/add-new-column/add-new-column.feature", import.meta.url);
  test("Add New Column on demog: the dialog, a formula built by hand, and the input history", {tag: ["@journey", "@realizes:powerpack.dialogs.add-new-column", "@realizes:powerpack.cp.add-new-column-persists"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await run.scenario("The toolbar icon opens the dialog, and its controls carry tooltips", async () => {
      await session.step(30, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(31, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(32, "And formula editor should be visible", () => shouldBe(page, el("formula editor"), "visible"));
      await session.step(33, "And formula hint should contain text \"Type '$' to select a column\"", () => shouldContainText(page, el("formula hint"), "Type '$' to select a column"));
      await session.step(34, "When user hovers over column name input", () => hoverOver(page, el("column name input")));
      await session.step(35, "Then tooltip should contain text \"Column name.\"", () => shouldContainText(page, el("tooltip"), "Column name."));
      await session.step(36, "When user hovers over column type input", () => hoverOver(page, el("column type input")));
      await session.step(37, "Then tooltip should contain text \"type is determined based on the expression\"", () => shouldContainText(page, el("tooltip"), "type is determined based on the expression"));
      await session.step(38, "When user hovers over preview grid viewer", () => hoverOver(page, el("preview grid viewer")));
      await session.step(39, "Then tooltip should contain text \"Preview result columns.\"", () => shouldContainText(page, el("tooltip"), "Preview result columns."));
      await session.step(40, "When user hovers over the \"cell 6 of __name\" area of column list viewer", () => hoverArea(page, "cell 6 of __name", el("column list viewer")));
      await session.step(41, "Then tooltip should contain text \"HEIGHT\"", () => shouldContainText(page, el("tooltip"), "HEIGHT"));
      await session.step(42, "When user hovers over functions sort icon", () => hoverOver(page, el("functions sort icon")));
      await session.step(43, "Then tooltip should contain text \"Select functions sort type\"", () => shouldContainText(page, el("tooltip"), "Select functions sort type"));
      await session.step(44, "When user hovers over name of \"Abs\" function entry", () => hoverOver(page, el("name of \"Abs\" function entry")));
      await session.step(45, "Then tooltip should contain text \"Abs\"", () => shouldContainText(page, el("tooltip"), "Abs"));
      await session.step(46, "When user hovers over \"History\" icon in \"Add New Column\" dialog", () => hoverOver(page, el("\"History\" icon in \"Add New Column\" dialog")));
      await session.step(47, "Then tooltip should contain text \"History\"", () => shouldContainText(page, el("tooltip"), "History"));
      await session.step(48, "When user hovers over \"Help\" icon in \"Add New Column\" dialog", () => hoverOver(page, el("\"Help\" icon in \"Add New Column\" dialog")));
      await session.step(49, "Then tooltip should contain text \"Help\"", () => shouldContainText(page, el("tooltip"), "Help"));
      await session.step(50, "When user hovers over \"Close\" icon in \"Add New Column\" dialog", () => hoverOver(page, el("\"Close\" icon in \"Add New Column\" dialog")));
      await session.step(51, "Then tooltip should contain text \"Close\"", () => shouldContainText(page, el("tooltip"), "Close"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The dialog resizes both ways with its content inside", async () => {
      await session.step(55, "When user remembers the size of \"Add New Column\" dialog", () => rememberSize(page, el("\"Add New Column\" dialog")));
      await session.step(56, "And user remembers the size of formula editor", () => rememberSize(page, el("formula editor")));
      await session.step(57, "And user remembers the size of preview grid viewer", () => rememberSize(page, el("preview grid viewer")));
      await session.step(58, "And user remembers the size of column list viewer", () => rememberSize(page, el("column list viewer")));
      await session.step(59, "And user remembers the size of functions panel", () => rememberSize(page, el("functions panel")));
      await session.step(60, "And user drags resize corner of Add New Column dialog by 200 and 150 pixels", () => dragCornerBy(page, el("resize corner of Add New Column dialog"), 200, 150));
      await session.step(61, "Then \"Add New Column\" dialog should be larger than remembered", () => sizeAgainstRemembered(page, el("\"Add New Column\" dialog"), "larger"));
      await session.step(62, "And formula editor should be wider than remembered", () => sizeAgainstRemembered(page, el("formula editor"), "wider"));
      await session.step(63, "And preview grid viewer should be larger than remembered", () => sizeAgainstRemembered(page, el("preview grid viewer"), "larger"));
      await session.step(64, "And column list viewer should keep its remembered width", () => keepsWidth(page, el("column list viewer")));
      await session.step(65, "And functions panel should keep its remembered width", () => keepsWidth(page, el("functions panel")));
      await session.step(66, "And formula editor should lie within \"Add New Column\" dialog", () => liesWithin(page, el("formula editor"), el("\"Add New Column\" dialog")));
      await session.step(67, "And column list viewer should lie within \"Add New Column\" dialog", () => liesWithin(page, el("column list viewer"), el("\"Add New Column\" dialog")));
      await session.step(68, "And functions panel should lie within \"Add New Column\" dialog", () => liesWithin(page, el("functions panel"), el("\"Add New Column\" dialog")));
      await session.step(69, "And preview grid viewer should lie within \"Add New Column\" dialog", () => liesWithin(page, el("preview grid viewer"), el("\"Add New Column\" dialog")));
      await session.step(70, "When user remembers the size of \"Add New Column\" dialog", () => rememberSize(page, el("\"Add New Column\" dialog")));
      await session.step(71, "And user remembers the size of formula editor", () => rememberSize(page, el("formula editor")));
      await session.step(72, "And user drags resize corner of Add New Column dialog by -300 and -200 pixels", () => dragCornerBy(page, el("resize corner of Add New Column dialog"), -300, -200));
      await session.step(73, "Then \"Add New Column\" dialog should be smaller than remembered", () => sizeAgainstRemembered(page, el("\"Add New Column\" dialog"), "smaller"));
      await session.step(74, "And formula editor should be narrower than remembered", () => sizeAgainstRemembered(page, el("formula editor"), "narrower"));
      await session.step(75, "And formula editor should lie within \"Add New Column\" dialog", () => liesWithin(page, el("formula editor"), el("\"Add New Column\" dialog")));
      await session.step(76, "And column list viewer should lie within \"Add New Column\" dialog", () => liesWithin(page, el("column list viewer"), el("\"Add New Column\" dialog")));
      await session.step(77, "And preview grid viewer should lie within \"Add New Column\" dialog", () => liesWithin(page, el("preview grid viewer"), el("\"Add New Column\" dialog")));
      await session.step(78, "And functions panel should lie within \"Add New Column\" dialog", () => liesWithin(page, el("functions panel"), el("\"Add New Column\" dialog")));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A formula from autocomplete and two dragged columns adds a column", async () => {
      await session.step(82, "When user types \"New\" into column name input", () => typeInto(page, "New", el("column name input")));
      await session.step(83, "And user types \"Rou\" into formula editor", () => typeInto(page, "Rou", el("formula editor")));
      await session.step(84, "Then completion list should be visible", () => shouldBe(page, el("completion list"), "visible"));
      await session.step(85, "When user accepts the highlighted completion with Enter", () => acceptCompletion(page, "Enter"));
      await session.step(86, "Then formula editor should hold the formula \"Round(a)\"", () => holdsFormula(page, el("formula editor"), "Round(a)"));
      await session.step(87, "And \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(88, "When user presses Delete in formula editor", () => pressKeyIn(page, "Delete", el("formula editor")));
      await session.step(89, "And user drags the \"header HEIGHT\" area of grid onto formula editor", () => dragAreaOnto(page, "header HEIGHT", el("grid"), el("formula editor")));
      await session.step(90, "Then formula editor should hold the formula \"Round(${HEIGHT})\"", () => holdsFormula(page, el("formula editor"), "Round(${HEIGHT})"));
      await session.step(91, "When user presses End in formula editor", () => pressKeyIn(page, "End", el("formula editor")));
      await session.step(92, "And user presses ArrowLeft in formula editor", () => pressKeyIn(page, "ArrowLeft", el("formula editor")));
      await session.step(93, "And user types \" + \" at the caret", () => typeAtCaret(page, " + "));
      await session.step(94, "And user drags the \"WEIGHT\" column of column list viewer onto formula editor", () => dragColumnName(page, "WEIGHT", el("column list viewer"), el("formula editor")));
      await session.step(95, "Then formula editor should hold the formula \"Round(${HEIGHT} + ${WEIGHT})\"", () => holdsFormula(page, el("formula editor"), "Round(${HEIGHT} + ${WEIGHT})"));
      await session.step(96, "And the preview grid should show \"New\" computed as numbers", () => previewComputes(page, "New"));
      await session.step(97, "When user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(98, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(99, "And the table should have a column \"New\"", () => hasColumn(page, "New"));
      await session.step(100, "And \"New\" column should have tag \"formula\" equal to \"Round(${HEIGHT} + ${WEIGHT})\"", () => columnTag(page, "New", "formula", "Round(${HEIGHT} + ${WEIGHT})"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The history of a reopened dialog fills the form back in", async () => {
      await session.step(104, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(105, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(106, "And column name input should have value \"\"", () => shouldHaveValue(page, el("column name input"), ""));
      await session.step(107, "And formula editor should hold the formula \"\"", () => holdsFormula(page, el("formula editor"), ""));
      await session.step(108, "When user clicks on \"History\" icon in \"Add New Column\" dialog", () => clickOn(page, el("\"History\" icon in \"Add New Column\" dialog")));
      await session.step(109, "Then input history menu should contain text \"Name: New\"", () => shouldContainText(page, el("input history menu"), "Name: New"));
      await session.step(110, "When user clicks on first menu item in input history menu", () => clickOn(page, el("first menu item in input history menu")));
      await session.step(111, "Then column name input should have value \"New\"", () => shouldHaveValue(page, el("column name input"), "New"));
      await session.step(112, "And formula editor should hold the formula \"Round(${HEIGHT} + ${WEIGHT})\"", () => holdsFormula(page, el("formula editor"), "Round(${HEIGHT} + ${WEIGHT})"));
      await session.step(113, "When user clicks on CANCEL button in \"Add New Column\" dialog", () => clickOn(page, el("CANCEL button in \"Add New Column\" dialog")));
      await session.step(114, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
