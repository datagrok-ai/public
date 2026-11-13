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
import {check, clickOn, pasteInto, selectIn, shouldBe, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnsExactly} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {balloonsClosed, contextPanelClosed, openDataset, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · PolyTool enumeration", () => {
  const session = feature(test, "features/demos/19-polytool-enumeration.feature", import.meta.url);
  test("Enumerate a peptide library", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(16, "And user opens HELM_sample dataset", () => openDataset(page, ds("HELM_sample")));
    await session.step(17, "Then \"HELM\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "HELM", "Macromolecule"));
    await session.step(18, "Given the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the \"HELM\" column of the grid is 900 px wide", () => gridColumnWidth(page, "HELM", 900));
    await session.step(20, "And the rows of the grid are 60 px tall", () => gridRowHeight(page, 60));
    await session.step(21, "When user clicks on the \"cell 1 of HELM\" area of grid", () => clickArea(page, "cell 1 of HELM", el("grid")));
    await session.step(22, "And user picks \"Bio > PolyTool > Enumerate HELM...\" from the top menu", () => pickFromTopMenu(page, "Bio > PolyTool > Enumerate HELM..."));
    await session.step(23, "Then \"PolyTool Helm Enumeration\" dialog should be visible", () => shouldBe(page, el("\"PolyTool Helm Enumeration\" dialog"), "visible"));
    await session.step(24, "When the video says \"Click the monomers to vary, then choose what to try there\"", () => videoSays(page, "Click the monomers to vary, then choose what to try there"));
    await session.step(25, "And user clicks on monomer 4 of the HELM in the \"PolyTool Helm Enumeration\" dialog", () => clickHelmMonomer(page, 4, "PolyTool Helm Enumeration"));
    await session.step(26, "And user double-clicks on the \"cell 1 of Monomers\" area of grid in Placeholders input in \"PolyTool Helm Enumeration\" dialog", () => doubleClickArea(page, "cell 1 of Monomers", el("grid in Placeholders input in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(27, "Then \"Select Monomers\" dialog should be visible", () => shouldBe(page, el("\"Select Monomers\" dialog"), "visible"));
    await session.step(28, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(29, "When user pastes \"D, Q, S\" into Search input in \"Select Monomers\" dialog", () => pasteInto(page, "D, Q, S", el("Search input in \"Select Monomers\" dialog")));
    await session.step(30, "And user clicks on OK button in \"Select Monomers\" dialog", () => clickOn(page, el("OK button in \"Select Monomers\" dialog")));
    await session.step(31, "And user clicks on monomer 12 of the HELM in the \"PolyTool Helm Enumeration\" dialog", () => clickHelmMonomer(page, 12, "PolyTool Helm Enumeration"));
    await session.step(32, "And user double-clicks on the \"cell 2 of Monomers\" area of grid in Placeholders input in \"PolyTool Helm Enumeration\" dialog", () => doubleClickArea(page, "cell 2 of Monomers", el("grid in Placeholders input in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(33, "Then \"Select Monomers\" dialog should be visible", () => shouldBe(page, el("\"Select Monomers\" dialog"), "visible"));
    await session.step(34, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(35, "When user pastes \"K, R\" into Search input in \"Select Monomers\" dialog", () => pasteInto(page, "K, R", el("Search input in \"Select Monomers\" dialog")));
    await session.step(36, "And user clicks on OK button in \"Select Monomers\" dialog", () => clickOn(page, el("OK button in \"Select Monomers\" dialog")));
    await session.step(37, "And user selects \"matrix\" in \"Enumerator type\" input in \"PolyTool Helm Enumeration\" dialog", () => selectIn(page, "matrix", el("\"Enumerator type\" input in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(38, "Then \"PolyTool Helm Enumeration\" dialog should contain the text \"6 sequences will be generated\"", () => shouldContainText(page, el("\"PolyTool Helm Enumeration\" dialog"), "6 sequences will be generated"));
    await session.step(39, "When user checks \"To atomic level\" checkbox in \"PolyTool Helm Enumeration\" dialog", () => check(page, el("\"To atomic level\" checkbox in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(40, "And user clicks on OK button in \"PolyTool Helm Enumeration\" dialog", () => clickOn(page, el("OK button in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(41, "Then the table should have 6 rows", () => rowCount(page, 6));
    await session.step(42, "And \"Enumerated\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "Enumerated", "Macromolecule"));
    await session.step(43, "And the table should have the columns \"Enumerated, molfile(Enumerated)\"", () => columnsExactly(page, "Enumerated, molfile(Enumerated)"));
    await session.step(44, "And the \"Table\" view should be current", () => viewIsCurrent(page, "Table"));
    await session.step(45, "Given the \"Enumerated\" column of the grid is 760 px wide", () => gridColumnWidth(page, "Enumerated", 760));
    await session.step(46, "And the \"molfile(Enumerated)\" column of the grid is 680 px wide", () => gridColumnWidth(page, "molfile(Enumerated)", 680));
    await session.step(47, "And the rows of the grid are 190 px tall", () => gridRowHeight(page, 190));
    await session.step(48, "When the video says \"Six peptides, each with its structure\"", () => videoSays(page, "Six peptides, each with its structure"));
    await session.step(49, "And user hovers over monomer 4 of the \"cell 3 of Enumerated\" area of grid", () => hoverMonomer(page, 4, "cell 3 of Enumerated", el("grid")));
    await session.step(50, "And the video says \"The monomer under the pointer lights up in the molecule\"", () => videoSays(page, "The monomer under the pointer lights up in the molecule"));
    await session.step(51, "Then no errors should have been logged", () => noErrors(page));
  });
});
