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
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And the Crux sketcher is open on \"C[C@H](O)[C@@H](C)C(=O)O\", showing R, S, E and Z labels", () => cruxOpenLabelled(page, "C[C@H](O)[C@@H](C)C(=O)O"));
    await session.step(16, "When user presses Control+A in Crux canvas", () => pressKeyIn(page, "Control+A", el("Crux canvas")), undefined, "Select everything (Ctrl+A)");
    await session.step(18, "And user opens the Crux context menu on the \"atom 1\" area", () => cruxMenu(page, "atom 1"), undefined, "Right-click a stereocentre");
    await session.step(20, "And user clicks on Crux enhanced stereo item", () => clickOn(page, el("Crux enhanced stereo item")), undefined, "Choose Enhanced Stereochemistry…");
    await session.step(22, "And user clicks on Crux AND button", () => clickOn(page, el("Crux AND button")), undefined, "Choose AND: the drawing and its mirror image, a racemic mixture");
    await session.step(24, "And user clicks on Crux enhanced stereo Apply button", () => clickOn(page, el("Crux enhanced stereo Apply button")), undefined, "Apply");
    await session.step(26, "Then the \"smiles\" reading of Crux sketcher widget should be \"C[C@H](O)[C@@H](C)C(=O)O |&1:1,3|\"", () => readingReads(page, "smiles", el("Crux sketcher widget"), "C[C@H](O)[C@@H](C)C(=O)O |&1:1,3|"), undefined, "Both centres are marked &1");
    await session.step(28, "When user opens the Crux context menu on the \"atom 1\" area", () => cruxMenu(page, "atom 1"), undefined, "Right-click a stereocentre again");
    await session.step(30, "And user clicks on Crux enhanced stereo item", () => clickOn(page, el("Crux enhanced stereo item")), undefined, "Choose Enhanced Stereochemistry…");
    await session.step(32, "And user clicks on Crux OR button", () => clickOn(page, el("Crux OR button")), undefined, "Choose OR: one of the two, not known which");
    await session.step(34, "And user clicks on Crux enhanced stereo Apply button", () => clickOn(page, el("Crux enhanced stereo Apply button")), undefined, "Apply");
    await session.step(36, "Then the \"smiles\" reading of Crux sketcher widget should be \"C[C@H](O)[C@@H](C)C(=O)O |o1:1,3|\"", () => readingReads(page, "smiles", el("Crux sketcher widget"), "C[C@H](O)[C@@H](C)C(=O)O |o1:1,3|"), undefined, "Both centres are marked or1");
  });
});
