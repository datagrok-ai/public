/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/19-polytool-enumeration.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {clickHelmMonomer} from '../../bindings/demo-polytool.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, pressKey, selectIn, shouldBe, shouldContainText, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnsExactly, setGridColumnFormat} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {balloonsClosed, contextPanelClosed, openDataset, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("5 · Design a library · Enumeration", () => {
  const session = feature(test, "features/demos/19-polytool-enumeration.feature", import.meta.url);
  test("Enumerate a peptide library", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(18, "And user opens HELM_sample dataset", () => openDataset(page, ds("HELM_sample")));
    await session.step(19, "Then \"HELM\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "HELM", "Macromolecule"));
    await session.step(20, "Given the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And user sets the format of \"Activity\" column in the grid to \"0.00\"", () => setGridColumnFormat(page, "Activity", "0.00"));
    await session.step(22, "And the \"HELM\" column of the grid is 900 px wide", () => gridColumnWidth(page, "HELM", 900));
    await session.step(23, "And the rows of the grid are 60 px tall", () => gridRowHeight(page, 60));
    await session.step(24, "When user clicks on the \"cell 1 of HELM\" area of grid", () => clickArea(page, "cell 1 of HELM", el("grid")));
    await session.step(25, "And user picks \"Bio > PolyTool > Enumerate HELM...\" from the top menu", () => pickFromTopMenu(page, "Bio > PolyTool > Enumerate HELM..."));
    await session.step(26, "Then \"PolyTool Helm Enumeration\" dialog should be visible", () => shouldBe(page, el("\"PolyTool Helm Enumeration\" dialog"), "visible"));
    await session.step(27, "When the video says \"Click the monomers to vary, then pick what to try\"", () => videoSays(page, "Click the monomers to vary, then pick what to try"));
    await session.step(28, "And user clicks on monomer 4 of the HELM in the \"PolyTool Helm Enumeration\" dialog", () => clickHelmMonomer(page, 4, "PolyTool Helm Enumeration"));
    await session.step(29, "And user double-clicks on the \"cell 1 of Monomers\" area of grid in Placeholders input in \"PolyTool Helm Enumeration\" dialog", () => doubleClickArea(page, "cell 1 of Monomers", el("grid in Placeholders input in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(30, "Then \"Select Monomers\" dialog should be visible", () => shouldBe(page, el("\"Select Monomers\" dialog"), "visible"));
    await session.step(31, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(32, "When user types \"D\" into Search input in \"Select Monomers\" dialog", () => typeInto(page, "D", el("Search input in \"Select Monomers\" dialog")));
    await session.step(33, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(34, "And user types \"Q\" into Search input in \"Select Monomers\" dialog", () => typeInto(page, "Q", el("Search input in \"Select Monomers\" dialog")));
    await session.step(35, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(36, "And user types \"S\" into Search input in \"Select Monomers\" dialog", () => typeInto(page, "S", el("Search input in \"Select Monomers\" dialog")));
    await session.step(37, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(38, "Then \"Select Monomers\" dialog should contain the text \"3 monomer(s) selected\"", () => shouldContainText(page, el("\"Select Monomers\" dialog"), "3 monomer(s) selected"));
    await session.step(39, "When user clicks on OK button in \"Select Monomers\" dialog", () => clickOn(page, el("OK button in \"Select Monomers\" dialog")));
    await session.step(40, "And user clicks on monomer 12 of the HELM in the \"PolyTool Helm Enumeration\" dialog", () => clickHelmMonomer(page, 12, "PolyTool Helm Enumeration"));
    await session.step(41, "And user double-clicks on the \"cell 2 of Monomers\" area of grid in Placeholders input in \"PolyTool Helm Enumeration\" dialog", () => doubleClickArea(page, "cell 2 of Monomers", el("grid in Placeholders input in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(42, "Then \"Select Monomers\" dialog should be visible", () => shouldBe(page, el("\"Select Monomers\" dialog"), "visible"));
    await session.step(43, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(44, "When user types \"K\" into Search input in \"Select Monomers\" dialog", () => typeInto(page, "K", el("Search input in \"Select Monomers\" dialog")));
    await session.step(45, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(46, "And user types \"R\" into Search input in \"Select Monomers\" dialog", () => typeInto(page, "R", el("Search input in \"Select Monomers\" dialog")));
    await session.step(47, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(48, "Then \"Select Monomers\" dialog should contain the text \"2 monomer(s) selected\"", () => shouldContainText(page, el("\"Select Monomers\" dialog"), "2 monomer(s) selected"));
    await session.step(49, "When user clicks on OK button in \"Select Monomers\" dialog", () => clickOn(page, el("OK button in \"Select Monomers\" dialog")));
    await session.step(50, "And user selects \"matrix\" in \"Enumerator type\" input in \"PolyTool Helm Enumeration\" dialog", () => selectIn(page, "matrix", el("\"Enumerator type\" input in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(51, "Then \"PolyTool Helm Enumeration\" dialog should contain the text \"6 sequences will be generated\"", () => shouldContainText(page, el("\"PolyTool Helm Enumeration\" dialog"), "6 sequences will be generated"));
    await session.step(52, "When user checks \"To atomic level\" checkbox in \"PolyTool Helm Enumeration\" dialog", () => check(page, el("\"To atomic level\" checkbox in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(53, "And user clicks on OK button in \"PolyTool Helm Enumeration\" dialog", () => clickOn(page, el("OK button in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(54, "Then the table should have 6 rows", () => rowCount(page, 6));
    await session.step(55, "And \"Enumerated\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "Enumerated", "Macromolecule"));
    await session.step(56, "And the table should have the columns \"Enumerated, molfile(Enumerated)\"", () => columnsExactly(page, "Enumerated, molfile(Enumerated)"));
    await session.step(57, "And the \"Table\" view should be current", () => viewIsCurrent(page, "Table"));
    await session.step(58, "Given the \"Enumerated\" column of the grid is 420 px wide", () => gridColumnWidth(page, "Enumerated", 420));
    await session.step(59, "And the \"molfile(Enumerated)\" column of the grid is 1100 px wide", () => gridColumnWidth(page, "molfile(Enumerated)", 1100));
    await session.step(60, "And the rows of the grid are 340 px tall", () => gridRowHeight(page, 340));
    await session.step(61, "When the video says \"Six variants, six molecules\"", () => videoSays(page, "Six variants, six molecules"));
    await session.step(62, "And user hovers over monomer 4 of the \"cell 2 of Enumerated\" area of grid", () => hoverMonomer(page, 4, "cell 2 of Enumerated", el("grid")));
    await session.step(63, "And the video says \"The monomer under the pointer lights up in the molecule\"", () => videoSays(page, "The monomer under the pointer lights up in the molecule"));
    await session.step(64, "Then no errors should have been logged", () => noErrors(page));
  });
});
