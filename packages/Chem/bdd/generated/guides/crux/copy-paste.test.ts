/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/copy-paste.feature
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
import {sketcherHolds} from '../../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, clipboardContains, pressKeyIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Copy and paste in Crux", () => {
  const session = feature(test, "features/guides/crux/copy-paste.feature", import.meta.url);
  test("Copy aspirin as SMILES, clear the canvas and paste it back", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And the Crux sketcher is open on \"CC(=O)Oc1ccccc1C(=O)O\"", () => cruxOpenOn(page, "CC(=O)Oc1ccccc1C(=O)O"));
    await session.step(16, "When user clicks on Crux copy as button", () => clickOn(page, el("Crux copy as button")), undefined, "Open Copy As");
    await session.step(18, "And user clicks on Crux SMILES item", () => clickOn(page, el("Crux SMILES item")), undefined, "Choose SMILES");
    await session.step(20, "Then the clipboard should contain text \"CC(=O)Oc1ccccc1C(=O)O\"", () => clipboardContains(page, "CC(=O)Oc1ccccc1C(=O)O"), undefined, "The clipboard holds aspirin's SMILES");
    await session.step(22, "When user clicks on Crux clear button", () => clickOn(page, el("Crux clear button")), undefined, "Clear the canvas");
    await session.step(24, "And user presses Control+V in Crux canvas", () => pressKeyIn(page, "Control+V", el("Crux canvas")), undefined, "Press Ctrl+V on the canvas");
    await session.step(26, "Then the sketcher in sketcher dialog should hold the molecule \"CC(=O)Oc1ccccc1C(=O)O\"", () => sketcherHolds(page, el("sketcher dialog"), "CC(=O)Oc1ccccc1C(=O)O"), undefined, "Aspirin is back");
  });
});
