/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/footguns/sketcher-warn.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '../../bindings/reaction-editor.js';
import '../../bindings/reaction-export.js';
import '../../bindings/reactions-app-text.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {engineGate} from '../../bindings/engine.js';
import {atomList, checksForCell, conformerTable, keeps3D, moleculeTable, noFlagSketcher, putChain, putMolecule, putQuery, readsMolecule2, sketcherButtonEnabled, slotShowsWord, slotShowsWords, slotShowsWordsNoFlag, v3000Of, valueContains} from '../../bindings/footguns.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Foot-guns in the molecule cell editor with Crux: what OK writes that breaks or changes meaning downstream", () => {
  const session = feature(test, "features/footguns/sketcher-warn.feature", import.meta.url);
  test("A query drawn into a molblock column is written, with a warn that calculators ignore it", {tag: ["@footguns", "@batch-2", "@catalogue-MOL-04"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(17, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(21, "Given user opens a table \"molecules\" with a molblock column \"molecule\" holding \"c1ccccc1\"", () => moleculeTable(page, "molecules", "molblock", "molecule", "c1ccccc1"));
    await session.step(22, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(23, "And user puts the query \"[C,N]1=CC=CC=C1\" into the Crux sketcher of the dialog", () => putQuery(page, "[C,N]1=CC=CC=C1"));
    await session.step(24, "Then the middle slot of the footer shows a warn flag mentioning \"query\"", () => slotShowsWord(page, "warn", "query"));
    await session.step(25, "And the OK button of the sketcher dialog is enabled", () => sketcherButtonEnabled(page, "OK"));
    await session.step(26, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(27, "Then the cell in row 1 of \"molecule\" holds a molblock with an atom list", () => atomList(page, 1, "molecule"));
    await session.step(28, "And the context panel's Checks section for that cell lists a warn mentioning \"query\"", () => checksForCell(page, "warn", "query"));
  });
  test("An edit of only enhanced stereo, an attachment point or an S-group is written, not thrown away [before=C[C@H](O)[C@@H](C)F, after=C[C@H](O)[C@@H](C)F |o1:1,3|, kept=o1:]", {tag: ["@footguns", "@batch-2", "@catalogue-MOL-07"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(17, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(32, "Given user opens a table \"molecules\" with a SMILES column \"molecule\" holding \"C[C@H](O)[C@@H](C)F\"", () => moleculeTable(page, "molecules", "SMILES", "molecule", "C[C@H](O)[C@@H](C)F"));
    await session.step(33, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(34, "And user puts \"C[C@H](O)[C@@H](C)F |o1:1,3|\" into the Crux sketcher of the dialog", () => putMolecule(page, "C[C@H](O)[C@@H](C)F |o1:1,3|"));
    await session.step(35, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(36, "Then the value in row 1 of \"molecule\" column contains \"o1:\"", () => valueContains(page, 1, "molecule", "o1:"));
    await session.step(37, "And the value in row 1 of \"molecule\" column is read by Chem as a molecule", () => readsMolecule2(page, 1, "molecule"));
  });
  test("Editing a 3D molfile in 2D warns that the z coordinates are lost; OK without an edit keeps them", {tag: ["@footguns", "@batch-2", "@catalogue-MOL-10"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(17, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(45, "Given user opens a table \"conformers\" with a molblock column \"molecule\" holding a 3D molfile of \"CCO\"", () => conformerTable(page, "conformers", "molecule", "CCO"));
    await session.step(46, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(47, "Then the footer of the sketcher dialog shows no flag", () => noFlagSketcher(page));
    await session.step(48, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(49, "Then the cell in row 1 of \"molecule\" still holds its 3D coordinates", () => keeps3D(page, 1, "molecule"));
    await session.step(50, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(51, "And user puts \"CCN\" into the Crux sketcher of the dialog", () => putMolecule(page, "CCN"));
    await session.step(52, "Then the middle slot of the footer shows a warn flag mentioning \"3D\"", () => slotShowsWord(page, "warn", "3D"));
    await session.step(53, "And the OK button of the sketcher dialog is enabled", () => sketcherButtonEnabled(page, "OK"));
  });
  test("Over 999 atoms in a V2000 column warns that the cell is written as V3000", {tag: ["@footguns", "@batch-2", "@catalogue-MOL-11"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(17, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(57, "Given user opens a table \"big\" with a molblock column \"molecule\" holding \"CCO\"", () => moleculeTable(page, "big", "molblock", "molecule", "CCO"));
    await session.step(58, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(59, "And user puts a chain of 1,200 carbon atoms into the Crux sketcher of the dialog", () => putChain(page, 1, 200));
    await session.step(60, "Then the middle slot of the footer shows a warn flag mentioning \"999\" and \"V3000\"", () => slotShowsWords(page, "warn", "999", "V3000"));
    await session.step(61, "And the OK button of the sketcher dialog is enabled", () => sketcherButtonEnabled(page, "OK"));
    await session.step(62, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(63, "Then the cell in row 1 of \"molecule\" holds a V3000 molblock of 1,200 atoms", () => v3000Of(page, 1, "molecule", 1, 200));
  });
  test("Two molecules drawn in one cell get a note, never a block", {tag: ["@footguns", "@batch-2", "@catalogue-MOL-13"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await session.step(17, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(67, "Given user opens a table \"molecules\" with a SMILES column \"molecule\" holding \"CCO\"", () => moleculeTable(page, "molecules", "SMILES", "molecule", "CCO"));
    await session.step(68, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(69, "And user puts \"CCO.c1ccccc1\" into the Crux sketcher of the dialog", () => putMolecule(page, "CCO.c1ccccc1"));
    await session.step(70, "Then the middle slot of the footer shows a note mentioning \"2\" and \"fragments\"", () => slotShowsWordsNoFlag(page, "note", "2", "fragments"));
    await session.step(71, "And the OK button of the sketcher dialog is enabled", () => sketcherButtonEnabled(page, "OK"));
  });
});
