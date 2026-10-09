/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/ez-double-bonds.feature
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
import {cruxMenu, cruxOpenLabelled} from '../../../bindings/crux.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("E and Z double bonds in Crux", () => {
  const session = feature(test, "features/guides/crux/ez-double-bonds.feature", import.meta.url);
  test("Make a double bond E, then Z, from its menu", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And the Crux sketcher is open on \"CC=CC(=O)O\", showing R, S, E and Z labels", () => cruxOpenLabelled(page, "CC=CC(=O)O"));
    await session.step(16, "When user opens the Crux context menu on the \"bond 1\" area", () => cruxMenu(page, "bond 1"), undefined, "Right-click the C=C double bond");
    await session.step(18, "And user clicks on Crux E item", () => clickOn(page, el("Crux E item")), undefined, "Choose E");
    await session.step(20, "Then the \"smiles\" reading of Crux sketcher widget should be the molecule \"C/C=C/C(=O)O\"", () => readingIsMolecule(page, "smiles", el("Crux sketcher widget"), "C/C=C/C(=O)O"), undefined, "(E)-crotonic acid: the methyl and the acid on opposite sides");
    await session.step(22, "When user opens the Crux context menu on the \"bond 1\" area", () => cruxMenu(page, "bond 1"), undefined, "Right-click the double bond again");
    await session.step(24, "And user clicks on Crux Z item", () => clickOn(page, el("Crux Z item")), undefined, "Choose Z");
    await session.step(26, "Then the \"smiles\" reading of Crux sketcher widget should be the molecule \"C/C=C\\C(=O)O\"", () => readingIsMolecule(page, "smiles", el("Crux sketcher widget"), "C/C=C\\C(=O)O"), undefined, "The bond flips to Z: isocrotonic acid");
  });
});
