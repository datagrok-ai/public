/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/footguns/sketcher-block.feature
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
import {moleculeTable, noFlagSketcher, putMolecule, putQuery, readsMolecule, sketcherButtonDisabled, sketcherButtonEnabled, slotShowsWord, slotShowsWords} from '../../bindings/footguns.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {rowMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {doubleClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Foot-guns in the molecule cell editor with Crux: what Chem cannot read back is blocked", () => {
  const session = feature(test, "features/footguns/sketcher-block.feature", import.meta.url);
  test("A valence Chem cannot read blocks OK in a molblock column [column kind=molblock, drawing=C(C)(C)(C)(C)C, keyword=valence]", {tag: ["@footguns", "@batch-1", "@catalogue-MOL-01"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(21, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(25, "Given user opens a table \"molecules\" with a molblock column \"molecule\" holding \"CCO\"", () => moleculeTable(page, "molecules", "molblock", "molecule", "CCO"));
    await session.step(26, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(27, "And user puts \"C(C)(C)(C)(C)C\" into the Crux sketcher of the dialog", () => putMolecule(page, "C(C)(C)(C)(C)C"));
    await session.step(28, "Then the middle slot of the footer shows a block flag mentioning \"valence\"", () => slotShowsWord(page, "block", "valence"));
    await session.step(29, "And the OK button of the sketcher dialog is disabled", () => sketcherButtonDisabled(page, "OK"));
    await session.step(30, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
    await session.step(31, "Then the molecule in row 1 of \"molecule\" column should be \"CCO\"", () => rowMolecule(page, 1, "molecule", "CCO"));
  });
  test("A valence Chem cannot read blocks OK in a SMILES column [column kind=SMILES, drawing=CN(C)(C)C, keyword=valence]", {tag: ["@footguns", "@batch-1", "@catalogue-MOL-01"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(21, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(25, "Given user opens a table \"molecules\" with a SMILES column \"molecule\" holding \"CCO\"", () => moleculeTable(page, "molecules", "SMILES", "molecule", "CCO"));
    await session.step(26, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(27, "And user puts \"CN(C)(C)C\" into the Crux sketcher of the dialog", () => putMolecule(page, "CN(C)(C)C"));
    await session.step(28, "Then the middle slot of the footer shows a block flag mentioning \"valence\"", () => slotShowsWord(page, "block", "valence"));
    await session.step(29, "And the OK button of the sketcher dialog is disabled", () => sketcherButtonDisabled(page, "OK"));
    await session.step(30, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
    await session.step(31, "Then the molecule in row 1 of \"molecule\" column should be \"CCO\"", () => rowMolecule(page, 1, "molecule", "CCO"));
  });
  test("Fixing the valence clears the block and OK writes a molecule every consumer reads", {tag: ["@footguns", "@batch-1", "@catalogue-MOL-01"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(21, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(40, "Given user opens a table \"molecules\" with a SMILES column \"molecule\" holding \"CCO\"", () => moleculeTable(page, "molecules", "SMILES", "molecule", "CCO"));
    await session.step(41, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(42, "And user puts \"C(C)(C)(C)(C)C\" into the Crux sketcher of the dialog", () => putMolecule(page, "C(C)(C)(C)(C)C"));
    await session.step(43, "Then the middle slot of the footer shows a block flag mentioning \"valence\"", () => slotShowsWord(page, "block", "valence"));
    await session.step(44, "When user puts \"CC(C)(C)C\" into the Crux sketcher of the dialog", () => putMolecule(page, "CC(C)(C)C"));
    await session.step(45, "Then the footer of the sketcher dialog shows no flag", () => noFlagSketcher(page));
    await session.step(46, "And the OK button of the sketcher dialog is enabled", () => sketcherButtonEnabled(page, "OK"));
    await session.step(47, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(48, "Then the molecule in row 1 of \"molecule\" column should be \"CC(C)(C)C\"", () => rowMolecule(page, 1, "molecule", "CC(C)(C)C"));
    await session.step(49, "And Chem reads row 1 of \"molecule\" as a molecule, not a query", () => readsMolecule(page, 1, "molecule"));
  });
  test("An aromatic ring with no Kekulé structure blocks OK [column kind=SMILES, drawing=c1ccnc1]", {tag: ["@footguns", "@batch-1", "@catalogue-MOL-02"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(21, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(53, "Given user opens a table \"molecules\" with a SMILES column \"molecule\" holding \"c1ccccc1\"", () => moleculeTable(page, "molecules", "SMILES", "molecule", "c1ccccc1"));
    await session.step(54, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(55, "And user puts \"c1ccnc1\" into the Crux sketcher of the dialog", () => putMolecule(page, "c1ccnc1"));
    await session.step(56, "Then the middle slot of the footer shows a block flag mentioning \"aromatic\" and \"Kekulé\"", () => slotShowsWords(page, "block", "aromatic", "Kekulé"));
    await session.step(57, "And the OK button of the sketcher dialog is disabled", () => sketcherButtonDisabled(page, "OK"));
    await session.step(58, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
    await session.step(59, "Then the molecule in row 1 of \"molecule\" column should be \"c1ccccc1\"", () => rowMolecule(page, 1, "molecule", "c1ccccc1"));
  });
  test("A query drawn into a SMILES column blocks OK instead of emptying the cell [query=[C,N]1=CC=CC=C1]", {tag: ["@footguns", "@batch-1", "@catalogue-MOL-03"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(21, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(67, "Given user opens a table \"molecules\" with a SMILES column \"molecule\" holding \"c1ccccc1\"", () => moleculeTable(page, "molecules", "SMILES", "molecule", "c1ccccc1"));
    await session.step(68, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(69, "And user puts the query \"[C,N]1=CC=CC=C1\" into the Crux sketcher of the dialog", () => putQuery(page, "[C,N]1=CC=CC=C1"));
    await session.step(70, "Then the middle slot of the footer shows a block flag mentioning \"query\" and \"SMILES\"", () => slotShowsWords(page, "block", "query", "SMILES"));
    await session.step(71, "And the OK button of the sketcher dialog is disabled", () => sketcherButtonDisabled(page, "OK"));
    await session.step(72, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
    await session.step(73, "Then the molecule in row 1 of \"molecule\" column should be \"c1ccccc1\"", () => rowMolecule(page, 1, "molecule", "c1ccccc1"));
  });
});
