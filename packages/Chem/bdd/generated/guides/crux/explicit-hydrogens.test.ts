/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/explicit-hydrogens.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/datasets.js';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {cruxOpenOn} from '../../../bindings/crux.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Explicit hydrogens in Crux", () => {
  const session = feature(test, "features/guides/crux/explicit-hydrogens.feature", import.meta.url);
  test("Draw L-alanine's hydrogens, then fold them back", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And the Crux sketcher is open on \"C[C@H](N)C(=O)O\"", () => cruxOpenOn(page, "C[C@H](N)C(=O)O"));
    await session.step(16, "When user clicks on Crux hydrogens button", () => clickOn(page, el("Crux hydrogens button")), undefined, "Open the Hydrogens menu");
    await session.step(18, "And user clicks on Crux add hydrogens item", () => clickOn(page, el("Crux add hydrogens item")), undefined, "Choose Add explicit hydrogens");
    await session.step(20, "Then the \"atoms\" reading of Crux sketcher widget should be 13", () => readingIs(page, "atoms", el("Crux sketcher widget"), 13), undefined, "All 13 atoms are drawn, the 7 hydrogens included");
    await session.step(21, "And the \"smiles\" reading of Crux sketcher widget should be the molecule \"C[C@H](N)C(=O)O\"", () => readingIsMolecule(page, "smiles", el("Crux sketcher widget"), "C[C@H](N)C(=O)O"));
    await session.step(23, "When user clicks on Crux hydrogens button", () => clickOn(page, el("Crux hydrogens button")), undefined, "Open the Hydrogens menu again");
    await session.step(25, "And user clicks on Crux remove hydrogens item", () => clickOn(page, el("Crux remove hydrogens item")), undefined, "Choose Remove explicit hydrogens");
    await session.step(27, "Then the \"atoms\" reading of Crux sketcher widget should be 6", () => readingIs(page, "atoms", el("Crux sketcher widget"), 6), undefined, "The hydrogens fold back into their atoms: still L-alanine");
    await session.step(28, "And the \"smiles\" reading of Crux sketcher widget should be the molecule \"C[C@H](N)C(=O)O\"", () => readingIsMolecule(page, "smiles", el("Crux sketcher widget"), "C[C@H](N)C(=O)O"));
  });
});
