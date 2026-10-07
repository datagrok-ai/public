/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sequence-translator/crux-polytool.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, browsePanelOpen, openTableOf, packageInstalled, sketcherIs, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes, readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Crux in PolyTool's molecule dialogs", () => {
  const session = feature(test, "features/sequence-translator/crux-polytool.feature", import.meta.url);
  test("A core drawn in Crux with R1 enables the Draw Core dialog's OK, and the dialog reports R1", {tag: ["@crux", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(17, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(22, "Given user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["CCO"]]), [["molecule"],["CCO"]]);
    await session.step(25, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(26, "When user picks \"Chem > Transform > Markush Enumeration...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Markush Enumeration..."));
    await session.step(27, "Then \"Markush Enumerator\" dialog should be visible", () => shouldBe(page, el("\"Markush Enumerator\" dialog"), "visible"));
    await session.step(28, "When user clicks on \"Cores: open sketcher\" button in \"Markush Enumerator\" dialog", () => clickOn(page, el("\"Cores: open sketcher\" button in \"Markush Enumerator\" dialog")));
    await session.step(29, "Then \"Draw Core\" dialog should be visible", () => shouldBe(page, el("\"Draw Core\" dialog"), "visible"));
    await session.step(30, "And OK button in \"Draw Core\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Draw Core\" dialog"), "disabled"));
    await session.step(31, "When user clicks on crux benzene tool", () => clickOn(page, el("crux benzene tool")));
    await session.step(32, "And user clicks on crux canvas", () => clickOn(page, el("crux canvas")));
    await session.step(33, "And user clicks on crux single bond tool", () => clickOn(page, el("crux single bond tool")));
    await session.step(34, "And user clicks on the \"atom 0\" area of crux sketcher widget", () => clickArea(page, "atom 0", el("crux sketcher widget")));
    await session.step(35, "And user clicks on crux R-group tool", () => clickOn(page, el("crux R-group tool")));
    await session.step(36, "And user clicks on the \"atom 6\" area of crux sketcher widget", () => clickArea(page, "atom 6", el("crux sketcher widget")));
    await session.step(37, "And user clicks on crux R1 button", () => clickOn(page, el("crux R1 button")));
    await session.step(38, "And user clicks on crux R-Group OK button", () => clickOn(page, el("crux R-Group OK button")));
    await session.step(39, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"[*:1]c1ccccc1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "[*:1]c1ccccc1"));
    await session.step(40, "And \"Draw Core\" dialog should contain text \"Detected R-groups: R1.\"", () => shouldContainText(page, el("\"Draw Core\" dialog"), "Detected R-groups: R1."));
    await session.step(41, "And OK button in \"Draw Core\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Draw Core\" dialog"), "enabled"));
    await session.step(42, "When user clicks on CANCEL button in \"Draw Core\" dialog", () => clickOn(page, el("CANCEL button in \"Draw Core\" dialog")));
    await session.step(43, "And user clicks on CANCEL button in \"Markush Enumerator\" dialog", () => clickOn(page, el("CANCEL button in \"Markush Enumerator\" dialog")));
  });
  test("Each of a reaction rule's three molecules, edited in Crux, keeps its R labels", {tag: ["@crux", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(17, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(47, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(48, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(49, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(50, "And Files---App-Data---SequenceTranslator tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---SequenceTranslator tree node inside browse tree")));
    await session.step(51, "And Files---App-Data---SequenceTranslator---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---SequenceTranslator---samples tree node inside browse tree")));
    await session.step(52, "When user double-clicks Files---App-Data---SequenceTranslator---samples---cyclized.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---SequenceTranslator---samples---cyclized.csv tree node inside browse tree")));
    await session.step(53, "Then the \"cyclized\" view should be current", () => viewIsCurrent(page, "cyclized"));
    await session.step(54, "And \"seqs\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "seqs", "Macromolecule"));
    await session.step(55, "When user picks \"Bio > PolyTool > Convert...\" from the top menu", () => pickFromTopMenu(page, "Bio > PolyTool > Convert..."));
    await session.step(56, "Then \"PolyTool Conversion\" dialog should be visible", () => shouldBe(page, el("\"PolyTool Conversion\" dialog"), "visible"));
    await session.step(57, "When user clicks on \"Edit rules\" icon in \"PolyTool Conversion\" dialog", () => clickOn(page, el("\"Edit rules\" icon in \"PolyTool Conversion\" dialog")));
    await session.step(58, "Then the \"Manage Polytool Rules - rules_example.json\" view should be current", () => viewIsCurrent(page, "Manage Polytool Rules - rules_example.json"));
    await session.step(59, "When user clicks on \"Reactions\" tab", () => clickOn(page, el("\"Reactions\" tab")));
    await session.step(60, "And user clicks on \"Add rule\" button", () => clickOn(page, el("\"Add rule\" button")));
    await session.step(61, "Then \"Add Reaction Rule\" dialog should be visible", () => shouldBe(page, el("\"Add Reaction Rule\" dialog"), "visible"));
    await session.step(62, "When user clicks on editor of \"First reactant\" input in \"Add Reaction Rule\" dialog", () => clickOn(page, el("editor of \"First reactant\" input in \"Add Reaction Rule\" dialog")));
    await session.step(63, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"[*:1]C\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "[*:1]C"));
    await session.step(64, "When user clicks on crux single bond tool", () => clickOn(page, el("crux single bond tool")));
    await session.step(65, "And user clicks on the \"atom 1\" area of crux sketcher widget", () => clickArea(page, "atom 1", el("crux sketcher widget")));
    await session.step(66, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"[*:1]CC\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "[*:1]CC"));
    await session.step(67, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(68, "And user clicks on editor of \"Second reactant\" input in \"Add Reaction Rule\" dialog", () => clickOn(page, el("editor of \"Second reactant\" input in \"Add Reaction Rule\" dialog")));
    await session.step(69, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"[*:2]C\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "[*:2]C"));
    await session.step(70, "When user clicks on crux single bond tool", () => clickOn(page, el("crux single bond tool")));
    await session.step(71, "And user clicks on the \"atom 1\" area of crux sketcher widget", () => clickArea(page, "atom 1", el("crux sketcher widget")));
    await session.step(72, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"[*:2]CC\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "[*:2]CC"));
    await session.step(73, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(74, "And user clicks on editor of Product input in \"Add Reaction Rule\" dialog", () => clickOn(page, el("editor of Product input in \"Add Reaction Rule\" dialog")));
    await session.step(75, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"[*:1]CC[*:2]\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "[*:1]CC[*:2]"));
    await session.step(76, "When user clicks on crux single bond tool", () => clickOn(page, el("crux single bond tool")));
    await session.step(77, "And user clicks on the \"atom 1\" area of crux sketcher widget", () => clickArea(page, "atom 1", el("crux sketcher widget")));
    await session.step(78, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"[*:1]C(C)C[*:2]\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "[*:1]C(C)C[*:2]"));
    await session.step(79, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(80, "And user clicks on editor of \"First reactant\" input in \"Add Reaction Rule\" dialog", () => clickOn(page, el("editor of \"First reactant\" input in \"Add Reaction Rule\" dialog")));
    await session.step(81, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"[*:1]CC\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "[*:1]CC"));
    await session.step(82, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
    await session.step(83, "And user clicks on editor of \"Second reactant\" input in \"Add Reaction Rule\" dialog", () => clickOn(page, el("editor of \"Second reactant\" input in \"Add Reaction Rule\" dialog")));
    await session.step(84, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"[*:2]CC\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "[*:2]CC"));
    await session.step(85, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
    await session.step(86, "And user clicks on editor of Product input in \"Add Reaction Rule\" dialog", () => clickOn(page, el("editor of Product input in \"Add Reaction Rule\" dialog")));
    await session.step(87, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"[*:1]C(C)C[*:2]\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "[*:1]C(C)C[*:2]"));
    await session.step(88, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
    await session.step(89, "And user clicks on CANCEL button in \"Add Reaction Rule\" dialog", () => clickOn(page, el("CANCEL button in \"Add Reaction Rule\" dialog")));
    await session.step(90, "Then \"Add Reaction Rule\" dialog should be absent", () => shouldBe(page, el("\"Add Reaction Rule\" dialog"), "absent"));
  });
});
