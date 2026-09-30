/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sequence-translator/oligo-polytool.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.menu.polytool.convert, bio.menu.polytool.combine-sequences]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnUnits, columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu, topMenuLists} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelOpen, dialogCloses, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {errorBalloonText, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bio | PolyTool on a custom-notation column: Convert and Combine Sequences", () => {
  const session = feature(test, "features/sequence-translator/oligo-polytool.feature", import.meta.url);
  test("Bio | PolyTool on a custom-notation column: Convert and Combine Sequences", {tag: ["@journey", "@realizes:bio.menu.polytool.convert", "@realizes:bio.menu.polytool.combine-sequences"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(17, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(18, "And Files---App-Data---SequenceTranslator tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---SequenceTranslator tree node inside browse tree")));
    await session.step(19, "And Files---App-Data---SequenceTranslator---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---SequenceTranslator---samples tree node inside browse tree")));
    await session.step(20, "When user double-clicks Files---App-Data---SequenceTranslator---samples---cyclized.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---SequenceTranslator---samples---cyclized.csv tree node inside browse tree")));
    await session.step(21, "Then the \"cyclized\" view should be current", () => viewIsCurrent(page, "cyclized"));
    await session.step(22, "And the table should have 14 rows", () => rowCount(page, 14));
    await session.step(23, "And \"seqs\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "seqs", "Macromolecule"));
    await run.scenario("The PolyTool submenu offers Convert, Enumerate HELM and Combine Sequences", async () => {
      await session.step(26, "Then the top menu should list:", () => topMenuLists(page, [["Bio > PolyTool > Convert..."],["Bio > PolyTool > Enumerate HELM..."],["Bio > PolyTool > Combine Sequences..."]]), [["Bio > PolyTool > Convert..."],["Bio > PolyTool > Enumerate HELM..."],["Bio > PolyTool > Combine Sequences..."]]);
    });
    await run.scenario("Convert with Get HELM adds a HELM column and a molfile column", async () => {
      await session.step(32, "When user picks \"Bio > PolyTool > Convert...\" from the top menu", () => pickFromTopMenu(page, "Bio > PolyTool > Convert..."));
      await session.step(33, "Then \"PolyTool Conversion\" dialog should be visible", () => shouldBe(page, el("\"PolyTool Conversion\" dialog"), "visible"));
      await session.step(34, "And Column input in \"PolyTool Conversion\" dialog should contain the text \"seqs\"", () => shouldContainText(page, el("Column input in \"PolyTool Conversion\" dialog"), "seqs"));
      await session.step(35, "And \"Get HELM\" checkbox in \"PolyTool Conversion\" dialog should be checked", () => shouldBe(page, el("\"Get HELM\" checkbox in \"PolyTool Conversion\" dialog"), "checked"));
      await session.step(36, "When user clicks on OK button in \"PolyTool Conversion\" dialog", () => clickOn(page, el("OK button in \"PolyTool Conversion\" dialog")));
      await session.step(37, "Then the \"PolyTool Conversion\" dialog should close", () => dialogCloses(page, "PolyTool Conversion"));
      await session.step(38, "And a new column \"transformed(seqs)\" should have been added", () => newColumnNamed(page, "transformed(seqs)"));
      await session.step(39, "And a new column \"molfile(seqs)\" should have been added", () => newColumnNamed(page, "molfile(seqs)"));
      await session.step(40, "And \"transformed(seqs)\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "transformed(seqs)", "Macromolecule"));
      await session.step(41, "And \"transformed(seqs)\" column should have units \"helm\"", () => columnUnits(page, "transformed(seqs)", "helm"));
      await session.step(42, "And \"molfile(seqs)\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molfile(seqs)", "Molecule"));
      await session.step(43, "And the value of \"transformed(seqs)\" column in row 1 should be \"PEPTIDE1{R.F.C.T.G.H.F.Y.G.H.F.Y.G.H.F.Y.P.C.[meI]}$PEPTIDE1,PEPTIDE1,3:R3-18:R3$$$V2.0\"", () => valueInRow(page, "transformed(seqs)", 1, "PEPTIDE1{R.F.C.T.G.H.F.Y.G.H.F.Y.G.H.F.Y.P.C.[meI]}$PEPTIDE1,PEPTIDE1,3:R3-18:R3$$$V2.0"));
      await session.step(44, "And \"transformed(seqs)\" column should have no missing values", () => columnComplete(page, "transformed(seqs)"));
      await session.step(45, "And \"molfile(seqs)\" column should have no missing values", () => columnComplete(page, "molfile(seqs)"));
      await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Combine Sequences with no table chosen says so and makes no table", async () => {
      await session.step(50, "When user picks \"Bio > PolyTool > Combine Sequences...\" from the top menu", () => pickFromTopMenu(page, "Bio > PolyTool > Combine Sequences..."));
      await session.step(51, "Then \"Combine Sequences\" dialog should be visible", () => shouldBe(page, el("\"Combine Sequences\" dialog"), "visible"));
      await session.step(52, "When user clicks on OK button in \"Combine Sequences\" dialog", () => clickOn(page, el("OK button in \"Combine Sequences\" dialog")));
      await session.step(53, "Then an error balloon containing \"Please fill all the fields\" should have been shown", () => errorBalloonText(page, "Please fill all the fields"));
      await session.step(54, "And the open table views should be exactly \"cyclized\"", () => openTableViewsExactly(page, "cyclized"));
    });
    await run.scenario("Combine Sequences of seqs with itself opens a table of all 196 pairs", async () => {
      await session.step(57, "When user picks \"Bio > PolyTool > Combine Sequences...\" from the top menu", () => pickFromTopMenu(page, "Bio > PolyTool > Combine Sequences..."));
      await session.step(58, "Then \"Combine Sequences\" dialog should be visible", () => shouldBe(page, el("\"Combine Sequences\" dialog"), "visible"));
      await session.step(59, "When user selects \"cyclized\" in Table input in \"Combine Sequences\" dialog", () => selectIn(page, "cyclized", el("Table input in \"Combine Sequences\" dialog")));
      await session.step(60, "Then Column input in \"Combine Sequences\" dialog should have the value \"seqs\"", () => shouldHaveValue(page, el("Column input in \"Combine Sequences\" dialog"), "seqs"));
      await session.step(61, "When user clicks on Add icon in \"Combine Sequences\" dialog", () => clickOn(page, el("Add icon in \"Combine Sequences\" dialog")));
      await session.step(62, "And user selects \"cyclized\" in second Table input in \"Combine Sequences\" dialog", () => selectIn(page, "cyclized", el("second Table input in \"Combine Sequences\" dialog")));
      await session.step(63, "Then second Column input in \"Combine Sequences\" dialog should have the value \"seqs\"", () => shouldHaveValue(page, el("second Column input in \"Combine Sequences\" dialog"), "seqs"));
      await session.step(64, "When user enters \"-\" into Separator input in \"Combine Sequences\" dialog", () => enterInto(page, "-", el("Separator input in \"Combine Sequences\" dialog")));
      await session.step(65, "And user clicks on OK button in \"Combine Sequences\" dialog", () => clickOn(page, el("OK button in \"Combine Sequences\" dialog")));
      await session.step(66, "Then the \"Combined Sequences\" view should be current", () => viewIsCurrent(page, "Combined Sequences"));
      await session.step(67, "And the table should have 196 rows", () => rowCount(page, 196));
      await session.step(68, "And the table should have the columns \"Combined Sequences\"", () => columnsExactly(page, "Combined Sequences"));
      await session.step(69, "And the value of \"Combined Sequences\" column in row 1 should be \"R-F-C(1)-T-G-H-F-Y-G-H-F-Y-G-H-F-Y-P-C(1)-meI-R-F-C(1)-T-G-H-F-Y-G-H-F-Y-G-H-F-Y-P-C(1)-meI\"", () => valueInRow(page, "Combined Sequences", 1, "R-F-C(1)-T-G-H-F-Y-G-H-F-Y-G-H-F-Y-P-C(1)-meI-R-F-C(1)-T-G-H-F-Y-G-H-F-Y-G-H-F-Y-P-C(1)-meI"));
      await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
