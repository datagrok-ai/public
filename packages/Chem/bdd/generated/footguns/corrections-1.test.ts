/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/footguns/corrections-1.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '../../bindings/reaction-editor.js';
import '../../bindings/reaction-export.js';
import '../../bindings/reactions-app-text.js';
import '../../bindings/sketcher-app.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {engineGate} from '../../bindings/engine.js';
import {clickDialogButton, clickReactionEditorButton, enterSmarts, eraseArrow, moleculeTable, noChecks, noFlagNamed, noFlagSketcher, noReactionPane, pageEngine, putBackReactions, putMolecule, reactionEditorButtonDisabled, reactionTable, savedTwoComponent, selectReaction, sketcherButtonEnabled, slotShowsWord, slotShowsWords, topButtonEnabled, valueAsOpened, valueChanged} from '../../bindings/footguns.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {setColumnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openTableOf, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea, doubleClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Corrections: RDKit mode on this page only, and Chem's own urea formation once fixed", () => {
  const session = feature(test, "features/footguns/corrections-1.feature", import.meta.url);
  test("RDKit mode shows no flag in the editor and lets OK write as before", {tag: ["@footguns", "@batch-6", "@catalogue-UI-03"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(19, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(24, "Given the Chemical Engine of this page is \"RDKit\"", () => pageEngine(page, "RDKit"));
    await session.step(25, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["CCO"]]), [["molecule"],["CCO"]]);
    await session.step(28, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(29, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(30, "And user puts \"C(C)(C)(C)(C)C\" into the Crux sketcher of the dialog", () => putMolecule(page, "C(C)(C)(C)(C)C"));
    await session.step(31, "Then the footer of the sketcher dialog shows no flag", () => noFlagSketcher(page));
    await session.step(32, "And the OK button of the sketcher dialog is enabled", () => sketcherButtonEnabled(page, "OK"));
    await session.step(33, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(34, "Then the value in row 1 of \"molecule\" column is no longer \"CCO\"", () => valueChanged(page, 1, "molecule", "CCO"));
  });
  test("RDKit mode: a template Crux mode flags has no Checks section and the reaction editor no flag", {tag: ["@footguns", "@batch-6", "@catalogue-UI-03"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(19, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(40, "Given the Chemical Engine of this page is \"RDKit\"", () => pageEngine(page, "RDKit"));
    await session.step(41, "And user opens a table \"templates\" with a reaction column \"template\" of units \"smarts\":", () => reactionTable(page, "templates", "template", "smarts", [["template"],["[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]"]]), [["template"],["[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]"]]);
    await session.step(44, "When user clicks on the \"cell 1 of template\" area of grid", () => clickArea(page, "cell 1 of template", el("grid")));
    await session.step(45, "Then the context panel has no Checks section", () => noChecks(page));
    await session.step(47, "Given user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["CCO"]]), [["molecule"],["CCO"]]);
    await session.step(50, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(51, "When user picks \"Chem | Transform | Reactions | Transformation...\" from the top menu", () => pickFromTopMenu(page, "Chem | Transform | Reactions | Transformation..."));
    await session.step(52, "And user clicks the \"New Reaction\" button of the dialog", () => clickDialogButton(page, "New Reaction"));
    await session.step(53, "And user enters \"[#6:1][C:2](=O)O>>[#6:1][C:2](=O)[*${n}]\" as the reaction SMARTS", () => enterSmarts(page, "[#6:1][C:2](=O)O>>[#6:1][C:2](=O)[*${n}]"));
    await session.step(54, "Then the footer of the \"New Reaction\" dialog shows no flag", () => noFlagNamed(page, "New Reaction"));
  });
  test("RDKit mode: the warns above show nowhere", {tag: ["@footguns", "@batch-6", "@catalogue-UI-03", "@catalogue-MOL-04", "@catalogue-MOL-05"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(19, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(59, "Given the Chemical Engine of this page is \"RDKit\"", () => pageEngine(page, "RDKit"));
    await session.step(60, "And user opens a table \"molecules\" with a SMILES column \"molecule\" holding \"CO\"", () => moleculeTable(page, "molecules", "SMILES", "molecule", "CO"));
    await session.step(61, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(62, "And user puts \"[CH3:1][OH:2]\" into the Crux sketcher of the dialog", () => putMolecule(page, "[CH3:1][OH:2]"));
    await session.step(63, "Then the footer of the sketcher dialog shows no flag", () => noFlagSketcher(page));
    await session.step(64, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(65, "And user clicks on the \"cell 1 of molecule\" area of grid", () => clickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(66, "Then the context panel has no Checks section", () => noChecks(page));
  });
  test("RDKit mode has no Reaction panel and no Checks section for these cells", {tag: ["@footguns", "@batch-6", "@catalogue-UI-03", "@catalogue-MAP-01"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(19, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(71, "Given the Chemical Engine of this page is \"RDKit\"", () => pageEngine(page, "RDKit"));
    await session.step(72, "And user opens a table \"reactions\" with a reaction column \"reaction\" of units \"smiles\":", () => reactionTable(page, "reactions", "reaction", "smiles", [["reaction"],["C(C)(C)(C)(C)C>>CO"]]), [["reaction"],["C(C)(C)(C)(C)C>>CO"]]);
    await session.step(75, "When user clicks on the \"cell 1 of reaction\" area of grid", () => clickArea(page, "cell 1 of reaction", el("grid")));
    await session.step(76, "Then the context panel has no Checks section", () => noChecks(page));
    await session.step(77, "And the context panel has no Reaction panel", () => noReactionPane(page));
  });
  test("RDKit mode: a double-click opens the reaction editor, and erasing the arrow blocks OK", {tag: ["@footguns", "@batch-6", "@catalogue-UI-03", "@catalogue-RXN-16"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(19, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(84, "Given the Chemical Engine of this page is \"RDKit\"", () => pageEngine(page, "RDKit"));
    await session.step(85, "And user opens a table \"reactions\" with a reaction column \"reaction\" of units \"smiles\":", () => reactionTable(page, "reactions", "reaction", "smiles", [["reaction"],["CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1"]]), [["reaction"],["CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1"]]);
    await session.step(88, "When user double-clicks on the \"cell 1 of reaction\" area of grid", () => doubleClickArea(page, "cell 1 of reaction", el("grid")));
    await session.step(89, "Then reaction editor dialog should be visible", () => shouldBe(page, el("reaction editor dialog"), "visible"));
    await session.step(90, "When user erases the arrow in the reaction editor", () => eraseArrow(page));
    await session.step(91, "Then the middle slot of the footer shows a block flag mentioning \"arrow\"", () => slotShowsWord(page, "block", "arrow"));
    await session.step(92, "And the OK button of the reaction editor is disabled", () => reactionEditorButtonDisabled(page, "OK"));
    await session.step(93, "When user clicks on CANCEL button of the reaction editor", () => clickReactionEditorButton(page, "CANCEL"));
    await session.step(94, "Then the value in row 1 of \"reaction\" column is exactly what it was", () => valueAsOpened(page, 1, "reaction"));
  });
  test("The Two-Component dialog warns before Run that a saved reaction takes three reactants", {tag: ["@footguns", "@batch-6", "@catalogue-TPL-10", "@catalogue-TRF-01"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(19, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(100, "Given user opens a table \"pairs\" with:", () => openTableOf(page, "pairs", [["acid","amine"],["CC(=O)O","CCN"]]), [["acid","amine"],["CC(=O)O","CCN"]]);
    await session.step(103, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(105, "And user sets the semantic type of \"acid\" column to \"Molecule\"", () => setColumnSemType(page, "acid", "Molecule"));
    await session.step(106, "And user sets the semantic type of \"amine\" column to \"Molecule\"", () => setColumnSemType(page, "amine", "Molecule"));
    await session.step(107, "And the user's saved reactions are put back when the feature ends", () => putBackReactions(page));
    await session.step(108, "And the user has a saved two-component reaction \"Three reactants\" \"[C:1](=[O:2])[OH].[NH2:3].CN(C)C=O>>[C:1](=[O:2])[NH:3]\"", () => savedTwoComponent(page, "Three reactants", "[C:1](=[O:2])[OH].[NH2:3].CN(C)C=O>>[C:1](=[O:2])[NH:3]"));
    await session.step(109, "When user picks \"Chem | Transform | Reactions | Two-Component Reaction...\" from the top menu", () => pickFromTopMenu(page, "Chem | Transform | Reactions | Two-Component Reaction..."));
    await session.step(110, "And user selects the reaction \"Three reactants\" in the dialog", () => selectReaction(page, "Three reactants"));
    await session.step(111, "Then the middle slot of the footer shows a warn flag mentioning \"3 reactants\" and \"2\"", () => slotShowsWords(page, "warn", "3 reactants", "2"));
    await session.step(112, "And the OK button of the dialog is enabled", () => topButtonEnabled(page, "OK"));
  });
});
