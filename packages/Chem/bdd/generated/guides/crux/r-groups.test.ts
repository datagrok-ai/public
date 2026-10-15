/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/r-groups.feature
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
import {check, clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("R-groups and attachment points in Crux", () => {
  const session = feature(test, "features/guides/crux/r-groups.feature", import.meta.url);
  test("Turn the phenol's OH into R1 and mark the acid's carbon as the attachment point", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And the Crux sketcher is open on \"Oc1ccc(cc1)C(=O)O\"", () => cruxOpenOn(page, "Oc1ccc(cc1)C(=O)O"));
    await session.step(15, "When user clicks on Crux R-group tool", () => clickOn(page, el("Crux R-group tool")), undefined, "Pick the R-group tool");
    await session.step(17, "And user clicks on the \"atom 0\" area of Crux sketcher widget", () => clickArea(page, "atom 0", el("Crux sketcher widget")), undefined, "Click the OH oxygen");
    await session.step(19, "And user clicks on Crux R1 button", () => clickOn(page, el("Crux R1 button")), undefined, "Choose R1");
    await session.step(21, "And user clicks on Crux R-Group OK button", () => clickOn(page, el("Crux R-Group OK button")), undefined, "OK");
    await session.step(23, "Then the \"smiles\" reading of Crux sketcher widget should be the molecule \"O=C(O)c1ccc([*:1])cc1\"", () => readingIsMolecule(page, "smiles", el("Crux sketcher widget"), "O=C(O)c1ccc([*:1])cc1"), undefined, "The oxygen is now R1");
    await session.step(25, "When user clicks on Crux R-group tool", () => clickOn(page, el("Crux R-group tool")), undefined, "Open the R-group palette");
    await session.step(27, "And user clicks on Crux attachment point tool", () => clickOn(page, el("Crux attachment point tool")), undefined, "Pick the attachment point tool");
    await session.step(29, "And user clicks on the \"atom 7\" area of Crux sketcher widget", () => clickArea(page, "atom 7", el("Crux sketcher widget")), undefined, "Click the acid's carbon");
    await session.step(31, "And user checks Crux primary attachment point checkbox", () => check(page, el("Crux primary attachment point checkbox")), undefined, "Mark it as the primary attachment point");
    await session.step(33, "And user clicks on Crux Attachment Points OK button", () => clickOn(page, el("Crux Attachment Points OK button")), undefined, "OK");
    await session.step(35, "Then the \"smiles\" reading of Crux sketcher widget should be \"O=C(O)c1ccc([*:1])cc1 |atomProp:1.molAttchpt.1:7.dummyLabel.R1|\"", () => readingReads(page, "smiles", el("Crux sketcher widget"), "O=C(O)c1ccc([*:1])cc1 |atomProp:1.molAttchpt.1:7.dummyLabel.R1|"), undefined, "The acid's carbon is marked as where the fragment attaches");
  });
});
