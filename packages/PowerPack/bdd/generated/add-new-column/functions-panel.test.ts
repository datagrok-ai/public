/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/functions-panel.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.dialogs.add-new-column, GROK-17004]
--- */
import {test} from '@playwright/test';
import '../../bindings/enrichment.js';
import '../../bindings/formula-lines.js';
import '../../bindings/home.js';
import '../../bindings/io.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clickColumnName, functionsAsRemembered, functionsByName, functionsNotAsRemembered, functionsStartWith, functionsTakeFirst, highlightStandsOut, highlightsExactly, highlightsInColor, holdsFormula, previewAbs, previewComputes, previewLacks, rememberFunctions} from '../../bindings/add-new-column.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, dragTo, hoverOver, pasteInto, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {menuLists, noErrors, pickFromOpenMenu, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The functions and columns of Add New Column on SPGI: insertion, auto-bound columns and sorting", () => {
  const session = feature(test, "features/add-new-column/functions-panel.feature", import.meta.url);
  test("The functions and columns of Add New Column on SPGI: insertion, auto-bound columns and sorting", {tag: ["@journey", "@realizes:powerpack.dialogs.add-new-column", "@realizes:GROK-17004"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(28, "And user opens SPGI dataset", () => openDataset(page, ds("SPGI")));
    await session.step(29, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
    await session.step(30, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
    await session.step(31, "And functions list should be visible", () => shouldBe(page, el("functions list"), "visible"));
    await run.scenario("With no column picked, the plus icon and a drag insert a function with its parameters", async () => {
      await session.step(34, "When user hovers over name of \"Abs\" function entry", () => hoverOver(page, el("name of \"Abs\" function entry")));
      await session.step(35, "Then plus of \"Abs\" function entry should be visible", () => shouldBe(page, el("plus of \"Abs\" function entry"), "visible"));
      await session.step(36, "When user clicks on plus of \"Abs\" function entry", () => clickOn(page, el("plus of \"Abs\" function entry")));
      await session.step(37, "Then formula editor should hold the formula \"Abs(x)\"", () => holdsFormula(page, el("formula editor"), "Abs(x)"));
      await session.step(38, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(39, "Then formula editor should hold the formula \"\"", () => holdsFormula(page, el("formula editor"), ""));
      await session.step(40, "And the preview grid should not show \"Abs(x)\"", () => previewLacks(page, "Abs(x)"));
      await session.step(41, "When user drags name of \"Abs\" function entry to formula editor", () => dragTo(page, el("name of \"Abs\" function entry"), el("formula editor")));
      await session.step(42, "Then formula editor should hold the formula \"Abs(x)\"", () => holdsFormula(page, el("formula editor"), "Abs(x)"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Structure brings the Molecule functions up, and getCLogP takes it", async () => {
      await session.step(46, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(47, "And user clicks on the \"Structure\" column in column list viewer", () => clickColumnName(page, "Structure", el("column list viewer")));
      await session.step(48, "Then the first 5 functions of the functions list should take a \"Molecule\" first", () => functionsTakeFirst(page, 5, "Molecule"));
      await session.step(49, "When user hovers over name of \"getCLogP\" function entry", () => hoverOver(page, el("name of \"getCLogP\" function entry")));
      await session.step(50, "And user clicks on plus of \"getCLogP\" function entry", () => clickOn(page, el("plus of \"getCLogP\" function entry")));
      await session.step(51, "Then formula editor should hold the formula \"Chem:getCLogP(${Structure})\"", () => holdsFormula(page, el("formula editor"), "Chem:getCLogP(${Structure})"));
      await session.step(52, "And the preview grid should show \"Chem:getCLogP(${Structure})\" computed as numbers", () => previewComputes(page, "Chem:getCLogP(${Structure})"));
      await session.step(53, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(54, "Then formula editor should hold the formula \"\"", () => holdsFormula(page, el("formula editor"), ""));
      await session.step(55, "And the preview grid should not show \"Chem:getCLogP(${Structure})\"", () => previewLacks(page, "Chem:getCLogP(${Structure})"));
      await session.step(56, "When user drags name of \"getCLogP\" function entry to formula editor", () => dragTo(page, el("name of \"getCLogP\" function entry"), el("formula editor")));
      await session.step(57, "Then formula editor should hold the formula \"Chem:getCLogP(${Structure})\"", () => holdsFormula(page, el("formula editor"), "Chem:getCLogP(${Structure})"));
      await session.step(58, "And the preview grid should show \"Chem:getCLogP(${Structure})\" computed as numbers", () => previewComputes(page, "Chem:getCLogP(${Structure})"));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A numeric column brings the numeric functions up, and Abs takes it", async () => {
      await session.step(62, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(63, "And user clicks on the \"Chemical Space X\" column in column list viewer", () => clickColumnName(page, "Chemical Space X", el("column list viewer")));
      await session.step(64, "Then the first 5 functions of the functions list should take a \"number\" first", () => functionsTakeFirst(page, 5, "number"));
      await session.step(65, "When user hovers over name of \"Abs\" function entry", () => hoverOver(page, el("name of \"Abs\" function entry")));
      await session.step(66, "And user clicks on plus of \"Abs\" function entry", () => clickOn(page, el("plus of \"Abs\" function entry")));
      await session.step(67, "Then formula editor should hold the formula \"Abs(${Chemical Space X})\"", () => holdsFormula(page, el("formula editor"), "Abs(${Chemical Space X})"));
      await session.step(68, "And the preview grid should show \"Abs(${Chemical Space X})\" as the absolute value of \"Chemical Space X\"", () => previewAbs(page, "Abs(${Chemical Space X})", "Chemical Space X"));
      await session.step(69, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(70, "Then formula editor should hold the formula \"\"", () => holdsFormula(page, el("formula editor"), ""));
      await session.step(71, "And the preview grid should not show \"Abs(${Chemical Space X})\"", () => previewLacks(page, "Abs(${Chemical Space X})"));
      await session.step(72, "When user drags name of \"Abs\" function entry to formula editor", () => dragTo(page, el("name of \"Abs\" function entry"), el("formula editor")));
      await session.step(73, "Then formula editor should hold the formula \"Abs(${Chemical Space X})\"", () => holdsFormula(page, el("formula editor"), "Abs(${Chemical Space X})"));
      await session.step(74, "And the preview grid should show \"Abs(${Chemical Space X})\" as the absolute value of \"Chemical Space X\"", () => previewAbs(page, "Abs(${Chemical Space X})", "Chemical Space X"));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A text column brings the text functions up", async () => {
      await session.step(78, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(79, "And user clicks on the \"Chemist\" column in column list viewer", () => clickColumnName(page, "Chemist", el("column list viewer")));
      await session.step(80, "Then the first 5 functions of the functions list should take a \"string\" first", () => functionsTakeFirst(page, 5, "string"));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("\"By name\" orders the functions alphabetically, and a column no longer reorders them", async () => {
      await session.step(84, "When user clicks on functions sort icon", () => clickOn(page, el("functions sort icon")));
      await session.step(85, "Then the open menu should list \"By name\"", () => menuLists(page, "By name"));
      await session.step(86, "And the open menu should list \"By relevance\"", () => menuLists(page, "By relevance"));
      await session.step(87, "When user picks \"By name\" from the open menu", () => pickFromOpenMenu(page, "By name"));
      await session.step(88, "Then the functions list should be sorted by name", () => functionsByName(page));
      await session.step(89, "And the functions list should start with \"Abs, Acos, Add\"", () => functionsStartWith(page, "Abs, Acos, Add"));
      await session.step(90, "When user remembers the order of the functions list", () => rememberFunctions(page));
      await session.step(91, "And user clicks on the \"Chemical Space X\" column in column list viewer", () => clickColumnName(page, "Chemical Space X", el("column list viewer")));
      await session.step(92, "Then the \"current row\" reading of column list viewer should be 19", () => readingIs(page, "current row", el("column list viewer"), 19));
      await session.step(93, "And the functions list should be in the remembered order", () => functionsAsRemembered(page));
      await session.step(94, "When user clicks on the \"Chemist\" column in column list viewer", () => clickColumnName(page, "Chemist", el("column list viewer")));
      await session.step(95, "Then the \"current row\" reading of column list viewer should be 5", () => readingIs(page, "current row", el("column list viewer"), 5));
      await session.step(96, "And the functions list should be in the remembered order", () => functionsAsRemembered(page));
      await session.step(97, "When user clicks on the \"Structure\" column in column list viewer", () => clickColumnName(page, "Structure", el("column list viewer")));
      await session.step(98, "Then the \"current row\" reading of column list viewer should be 2", () => readingIs(page, "current row", el("column list viewer"), 2));
      await session.step(99, "And the functions list should be in the remembered order", () => functionsAsRemembered(page));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column no parameter takes is not passed to the function", async () => {
      await session.step(103, "When user clicks on the \"Id\" column in column list viewer", () => clickColumnName(page, "Id", el("column list viewer")));
      await session.step(104, "Then the \"current row\" reading of column list viewer should be 1", () => readingIs(page, "current row", el("column list viewer"), 1));
      await session.step(105, "When user hovers over name of \"Abs\" function entry", () => hoverOver(page, el("name of \"Abs\" function entry")));
      await session.step(106, "And user clicks on plus of \"Abs\" function entry", () => clickOn(page, el("plus of \"Abs\" function entry")));
      await session.step(107, "Then formula editor should hold the formula \"Abs(x)\"", () => holdsFormula(page, el("formula editor"), "Abs(x)"));
      await session.step(108, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(109, "Then formula editor should hold the formula \"\"", () => holdsFormula(page, el("formula editor"), ""));
      await session.step(110, "And the preview grid should not show \"Abs(x)\"", () => previewLacks(page, "Abs(x)"));
      await session.step(111, "When user drags name of \"Abs\" function entry to formula editor", () => dragTo(page, el("name of \"Abs\" function entry"), el("formula editor")));
      await session.step(112, "Then formula editor should hold the formula \"Abs(x)\"", () => holdsFormula(page, el("formula editor"), "Abs(x)"));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("\"By relevance\" puts the column's functions back on top", async () => {
      await session.step(116, "When user clears formula editor", () => clearField(page, el("formula editor")));
      await session.step(117, "And user remembers the order of the functions list", () => rememberFunctions(page));
      await session.step(118, "And user clicks on functions sort icon", () => clickOn(page, el("functions sort icon")));
      await session.step(119, "And user picks \"By relevance\" from the open menu", () => pickFromOpenMenu(page, "By relevance"));
      await session.step(120, "And user clicks on the \"Structure\" column in column list viewer", () => clickColumnName(page, "Structure", el("column list viewer")));
      await session.step(121, "Then the functions list should not be in the remembered order", () => functionsNotAsRemembered(page));
      await session.step(122, "And the first 5 functions of the functions list should take a \"Molecule\" first", () => functionsTakeFirst(page, 5, "Molecule"));
      await session.step(123, "When user clicks on CANCEL button in \"Add New Column\" dialog", () => clickOn(page, el("CANCEL button in \"Add New Column\" dialog")));
      await session.step(124, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(125, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The formula of GROK-17004 pasted whole highlights all its column references", async () => {
      await session.step(128, "When user picks \"Edit > Add New Column...\" from the top menu", () => pickFromTopMenu(page, "Edit > Add New Column..."));
      await session.step(129, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(130, "When user pastes \"if(${Whole blood assay 1} != null, ${Whole blood assay 1}, if(${Route Admin}==\\\"PO\\\", ${Whole blood assay 1} / ${Chemical Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Species}, 'Rat') || Contains(${Species}, 'Rat Legacy'), 80, if(Contains(${Species}, 'Mouse'), 125, if(${Species}==\\\"Dog\\\", 30.9, if(${Species}==\\\"Monkey\\\", 43.6, if(${Species}==\\\"Minipig\\\", 39, null)))))*100\" into formula editor", () => pasteInto(page, "if(${Whole blood assay 1} != null, ${Whole blood assay 1}, if(${Route Admin}==\"PO\", ${Whole blood assay 1} / ${Chemical Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Species}, 'Rat') || Contains(${Species}, 'Rat Legacy'), 80, if(Contains(${Species}, 'Mouse'), 125, if(${Species}==\"Dog\", 30.9, if(${Species}==\"Monkey\", 43.6, if(${Species}==\"Minipig\", 39, null)))))*100", el("formula editor")));
      await session.step(131, "Then formula editor should contain text \"if(${Whole blood assay 1} != null\"", () => shouldContainText(page, el("formula editor"), "if(${Whole blood assay 1} != null"));
      await session.step(132, "And formula editor should contain text \"39, null)))))*100\"", () => shouldContainText(page, el("formula editor"), "39, null)))))*100"));
      await session.step(133, "And formula editor should highlight the column references \"${Whole blood assay 1}, ${Whole blood assay 1}, ${Route Admin}, ${Whole blood assay 1}, ${Chemical Space X}, ${Average Mass}, ${Species}, ${Species}, ${Species}, ${Species}, ${Species}, ${Species}\"", () => highlightsExactly(page, el("formula editor"), "${Whole blood assay 1}, ${Whole blood assay 1}, ${Route Admin}, ${Whole blood assay 1}, ${Chemical Space X}, ${Average Mass}, ${Species}, ${Species}, ${Species}, ${Species}, ${Species}, ${Species}"));
      await session.step(134, "And every column reference of formula editor should be drawn in the color of \"--blue-2\"", () => highlightsInColor(page, el("formula editor"), "--blue-2"));
      await session.step(135, "And every column reference of formula editor should differ in color from the plain text of its line", () => highlightStandsOut(page, el("formula editor")));
      await session.step(136, "And no errors should have been logged", () => noErrors(page));
      await session.step(137, "When user clicks on CANCEL button in \"Add New Column\" dialog", () => clickOn(page, el("CANCEL button in \"Add New Column\" dialog")));
      await session.step(138, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
