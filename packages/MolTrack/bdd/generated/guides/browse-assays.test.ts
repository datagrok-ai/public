/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/browse-assays.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {demoData, openEmsView} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {doubleClickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack: browse assays, their runs and results", () => {
  const session = feature(test, "features/guides/browse-assays.feature", import.meta.url);
  test("From the assay to a run, and to one of its results", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the MolTrack demo data is on the server", () => demoData(page));
    await session.step(12, "And user opens the \"Assays\" view of MolTrack", () => openEmsView(page, "Assays"));
    enter(page, "MolTrack");
    await session.step(14, "When user double-clicks on \"Hepatocyte Stability\" item in domain list", () => doubleClickOn(page, el("\"Hepatocyte Stability\" item in domain list")), undefined, "The assays");
    await session.step(16, "Then \"Assay runs\" tab should be visible", () => shouldBe(page, el("\"Assay runs\" tab"), "visible"), undefined, "The assay, and its runs below it");
    await session.step(17, "When user opens the \"Assay Results\" view of MolTrack", () => openEmsView(page, "Assay Results"));
    enter(page, "MolTrack");
    await session.step(19, "Then domain list should contain text \"uM\"", () => shouldContainText(page, el("domain list"), "uM"), undefined, "Every result: its batch, run, dose and CLint");
    await session.step(20, "When user double-clicks on first item in domain list", () => doubleClickOn(page, el("first item in domain list")));
    await session.step(22, "Then domain form should contain text \"Reported CLint\"", () => shouldContainText(page, el("domain form"), "Reported CLint"), undefined, "One result: the batch, the run, and the measured values");
  });
});
