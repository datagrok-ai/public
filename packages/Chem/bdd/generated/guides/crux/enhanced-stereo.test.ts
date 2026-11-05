/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/enhanced-stereo.feature
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
import {clickOn, pressKeyIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Enhanced stereochemistry in Crux", () => {
  const session = feature(test, "features/guides/crux/enhanced-stereo.feature", import.meta.url);
  test("Put two stereocentres in an AND group, then an OR group", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And the Crux sketcher is open on \"C[C@H](O)[C@@H](C)C(=O)O\", showing R, S, E and Z labels", () => cruxOpenLabelled(page, "C[C@H](O)[C@@H](C)C(=O)O"));
    await session.step(17, "When user presses Control+A in Crux canvas", () => pressKeyIn(page, "Control+A", el("Crux canvas")), undefined, "Select everything (Ctrl+A)");
    await session.step(19, "And user opens the Crux context menu on the \"atom 1\" area", () => cruxMenu(page, "atom 1"), undefined, "Right-click a stereocentre");
    await session.step(21, "And user clicks on Crux enhanced stereo item", () => clickOn(page, el("Crux enhanced stereo item")), undefined, "Choose Enhanced Stereochemistry…: the stereochemistry panel opens, both centres selected");
    await session.step(23, "And user clicks on Crux group chip of atom 1", () => clickOn(page, el("Crux group chip of atom 1")), undefined, "Open the group chip of a selected centre");
    await session.step(25, "And user clicks on Crux New AND group item", () => clickOn(page, el("Crux New AND group item")), undefined, "Choose a new AND group: the drawing and its mirror image, a racemic mixture");
    await session.step(27, "Then the \"smiles\" reading of Crux sketcher widget should be \"C[C@H](O)[C@@H](C)C(=O)O |&1:1,3|\"", () => readingReads(page, "smiles", el("Crux sketcher widget"), "C[C@H](O)[C@@H](C)C(=O)O |&1:1,3|"), undefined, "Both centres are marked &1");
    await session.step(29, "When user opens the Crux context menu on the \"atom 1\" area", () => cruxMenu(page, "atom 1"), undefined, "Right-click a stereocentre again");
    await session.step(31, "And user clicks on Crux enhanced stereo item", () => clickOn(page, el("Crux enhanced stereo item")), undefined, "Choose Enhanced Stereochemistry…: the panel, both centres still selected");
    await session.step(33, "And user clicks on Crux group chip of atom 1", () => clickOn(page, el("Crux group chip of atom 1")), undefined, "Open the group chip again");
    await session.step(35, "And user clicks on Crux New OR group item", () => clickOn(page, el("Crux New OR group item")), undefined, "Choose a new OR group: one of the two, not known which");
    await session.step(37, "Then the \"smiles\" reading of Crux sketcher widget should be \"C[C@H](O)[C@@H](C)C(=O)O |o1:1,3|\"", () => readingReads(page, "smiles", el("Crux sketcher widget"), "C[C@H](O)[C@@H](C)C(=O)O |o1:1,3|"), undefined, "Both centres are marked or1");
  });
});
