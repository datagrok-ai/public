/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/wedges-and-cip-labels.feature
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
import {cruxMenu, cruxOpenOn} from '../../../bindings/crux.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Wedges and R/S labels in Crux", () => {
  const session = feature(test, "features/guides/crux/wedges-and-cip-labels.feature", import.meta.url);
  test("Wedge a bond, read the centre's R or S, and invert it with a hashed wedge", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And the Crux sketcher is open on \"CC(N)C(=O)O\"", () => cruxOpenOn(page, "CC(N)C(=O)O"));
    await session.step(17, "When user clicks on Crux settings button", () => clickOn(page, el("Crux settings button")), undefined, "Open the gear menu");
    await session.step(19, "And user clicks on Crux settings item", () => clickOn(page, el("Crux settings item")), undefined, "Choose Settings…");
    await session.step(21, "And user checks Crux stereo labels checkbox", () => check(page, el("Crux stereo labels checkbox")), undefined, "Turn on \"Show R, S, E and Z labels\"");
    await session.step(23, "And user clicks on Crux settings Apply button", () => clickOn(page, el("Crux settings Apply button")), undefined, "Apply");
    await session.step(25, "And user clicks on Crux stereo bond tool", () => clickOn(page, el("Crux stereo bond tool")), undefined, "Pick the wedge bond tool");
    await session.step(27, "And user clicks on the \"bond 1\" area of Crux sketcher widget", () => clickArea(page, "bond 1", el("Crux sketcher widget")), undefined, "Click the bond from the centre to the nitrogen: it becomes a wedge");
    await session.step(29, "Then the \"smiles\" reading of Crux sketcher widget should be the molecule \"C[C@@H](N)C(=O)O\"", () => readingIsMolecule(page, "smiles", el("Crux sketcher widget"), "C[C@@H](N)C(=O)O"), undefined, "The centre is labelled (R): this is D-alanine");
    await session.step(31, "When user opens the Crux context menu on the \"bond 1\" area", () => cruxMenu(page, "bond 1"), undefined, "Right-click the wedge");
    await session.step(33, "And user clicks on Crux hash item", () => clickOn(page, el("Crux hash item")), undefined, "Choose Hashed wedge bond");
    await session.step(35, "Then the \"smiles\" reading of Crux sketcher widget should be the molecule \"C[C@H](N)C(=O)O\"", () => readingIsMolecule(page, "smiles", el("Crux sketcher widget"), "C[C@H](N)C(=O)O"), undefined, "The configuration flips from R to S: L-alanine");
  });
});
