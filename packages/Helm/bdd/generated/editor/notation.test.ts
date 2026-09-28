/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/editor/notation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [helm.cell-editor.molecule]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {helmInitialized, replaceNotation} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, shouldBe, shouldContainText, shouldHaveText, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColors, doubleClickArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Editing a sequence in the HELM editor", () => {
  const session = feature(test, "features/editor/notation.feature", import.meta.url);
  test("Editing a sequence in the HELM editor", {tag: ["@journey", "@realizes:helm.cell-editor.molecule"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the Helm package is initialized", () => helmInitialized(page));
    await session.step(18, "And user opens helm-showcase dataset", () => openDataset(page, ds("helm-showcase")));
    await session.step(19, "Then \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
    await run.scenario("The HELM tab shows the notation and an invalid edit shows a parse error", async () => {
      await session.step(22, "When user double-clicks on the \"cell 1 of HELM\" area of grid", () => doubleClickArea(page, "cell 1 of HELM", el("grid")));
      await session.step(23, "Then HELM editor should be visible", () => shouldBe(page, el("HELM editor"), "visible"));
      await session.step(24, "When user clicks on HELM tab", () => clickOn(page, el("HELM tab")));
      await session.step(25, "Then notation pane should have text \"PEPTIDE1{A.C}$$$$V2.0\"", () => shouldHaveText(page, el("notation pane"), "PEPTIDE1{A.C}$$$$V2.0"));
      await session.step(26, "And notation error should have text \"\"", () => shouldHaveText(page, el("notation error"), ""));
      await session.step(27, "When user replaces the text of notation pane with \"PEPTIDE1{A.ZZZ\"", () => replaceNotation(page, el("notation pane"), "PEPTIDE1{A.ZZZ"));
      await session.step(28, "And user presses Enter in notation pane", () => pressKeyIn(page, "Enter", el("notation pane")));
      await session.step(29, "Then notation error should contain text \"Expected '}' to close polymer 'PEPTIDE1'\"", () => shouldContainText(page, el("notation error"), "Expected '}' to close polymer 'PEPTIDE1'"));
      await session.step(30, "And there should be 2 visible drawn monomers", () => visibleCount(page, 2, el("drawn monomers")));
      await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A valid notation clears the error and redraws, without touching the cell", async () => {
      await session.step(34, "When user replaces the text of notation pane with \"PEPTIDE1{A.C.G}$$$$V2.0\"", () => replaceNotation(page, el("notation pane"), "PEPTIDE1{A.C.G}$$$$V2.0"));
      await session.step(35, "And user presses Enter in notation pane", () => pressKeyIn(page, "Enter", el("notation pane")));
      await session.step(36, "Then notation error should have text \"\"", () => shouldHaveText(page, el("notation error"), ""));
      await session.step(37, "And there should be 3 visible drawn monomers", () => visibleCount(page, 3, el("drawn monomers")));
      await session.step(38, "And the value of \"HELM\" column in row 1 should be \"PEPTIDE1{A.C}$$$$\"", () => valueInRow(page, "HELM", 1, "PEPTIDE1{A.C}$$$$"));
      await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Undo and redo step over the edit", async () => {
      await session.step(43, "When user clicks on undo button", () => clickOn(page, el("undo button")));
      await session.step(44, "Then there should be 2 visible drawn monomers", () => visibleCount(page, 2, el("drawn monomers")));
      await session.step(45, "And notation pane should have text \"PEPTIDE1{A.C}$$$$V2.0\"", () => shouldHaveText(page, el("notation pane"), "PEPTIDE1{A.C}$$$$V2.0"));
      await session.step(46, "When user clicks on redo button", () => clickOn(page, el("redo button")));
      await session.step(47, "Then there should be 3 visible drawn monomers", () => visibleCount(page, 3, el("drawn monomers")));
      await session.step(48, "And notation pane should have text \"PEPTIDE1{A.C.G}$$$$V2.0\"", () => shouldHaveText(page, el("notation pane"), "PEPTIDE1{A.C.G}$$$$V2.0"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clean layout keeps the structure", async () => {
      await session.step(52, "When user clicks on clean layout button", () => clickOn(page, el("clean layout button")));
      await session.step(53, "Then there should be 3 visible drawn monomers", () => visibleCount(page, 3, el("drawn monomers")));
      await session.step(54, "And notation pane should have text \"PEPTIDE1{A.C.G}$$$$V2.0\"", () => shouldHaveText(page, el("notation pane"), "PEPTIDE1{A.C.G}$$$$V2.0"));
      await session.step(55, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cancel discards the edit", async () => {
      await session.step(59, "When user clicks on CANCEL button in HELM editor", () => clickOn(page, el("CANCEL button in HELM editor")));
      await session.step(60, "Then HELM editor should be absent", () => shouldBe(page, el("HELM editor"), "absent"));
      await session.step(61, "And the value of \"HELM\" column in row 1 should be \"PEPTIDE1{A.C}$$$$\"", () => valueInRow(page, "HELM", 1, "PEPTIDE1{A.C}$$$$"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Properties tab computes the formula and the molecular weight", async () => {
      await session.step(65, "When user double-clicks on the \"cell 1 of HELM\" area of grid", () => doubleClickArea(page, "cell 1 of HELM", el("grid")));
      await session.step(66, "Then HELM editor should be visible", () => shouldBe(page, el("HELM editor"), "visible"));
      await session.step(67, "When user clicks on Properties tab", () => clickOn(page, el("Properties tab")));
      await session.step(68, "Then formula field should have text \"C6H12N2O3S\"", () => shouldHaveText(page, el("formula field"), "C6H12N2O3S"));
      await session.step(69, "And molecular weight field should have text \"192.23\"", () => shouldHaveText(page, el("molecular weight field"), "192.23"));
      await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("OK writes a sequence trimmed by its last monomer to the cell", async () => {
      await session.step(73, "When user clicks on CANCEL button in HELM editor", () => clickOn(page, el("CANCEL button in HELM editor")));
      await session.step(74, "And user double-clicks on the \"cell 2 of HELM\" area of grid", () => doubleClickArea(page, "cell 2 of HELM", el("grid")));
      await session.step(75, "Then HELM editor should be visible", () => shouldBe(page, el("HELM editor"), "visible"));
      await session.step(76, "And there should be 10 visible drawn monomers", () => visibleCount(page, 10, el("drawn monomers")));
      await session.step(77, "When user clicks on HELM tab", () => clickOn(page, el("HELM tab")));
      await session.step(78, "And user replaces the text of notation pane with \"PEPTIDE1{A.C.D.E.F.G.H.I.K}$$$$V2.0\"", () => replaceNotation(page, el("notation pane"), "PEPTIDE1{A.C.D.E.F.G.H.I.K}$$$$V2.0"));
      await session.step(79, "And user presses Enter in notation pane", () => pressKeyIn(page, "Enter", el("notation pane")));
      await session.step(80, "Then there should be 9 visible drawn monomers", () => visibleCount(page, 9, el("drawn monomers")));
      await session.step(81, "And the value of \"HELM\" column in row 2 should be \"PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$\"", () => valueInRow(page, "HELM", 2, "PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$"));
      await session.step(82, "When user clicks on OK button in HELM editor", () => clickOn(page, el("OK button in HELM editor")));
      await session.step(83, "Then HELM editor should be absent", () => shouldBe(page, el("HELM editor"), "absent"));
      await session.step(84, "And the value of \"HELM\" column in row 2 should be \"PEPTIDE1{A.C.D.E.F.G.H.I.K}$$$$V2.0\"", () => valueInRow(page, "HELM", 2, "PEPTIDE1{A.C.D.E.F.G.H.I.K}$$$$V2.0"));
      await session.step(85, "And the \"cell 2 of HELM\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 2 of HELM", el("grid"), 3));
      await session.step(86, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
