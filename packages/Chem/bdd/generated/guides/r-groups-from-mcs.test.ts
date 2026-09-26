/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/r-groups-from-mcs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, finishedUpdating, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openDataset, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Break a compound series into its core and R-groups", () => {
  const session = feature(test, "features/guides/r-groups-from-mcs.feature", import.meta.url);
  test("Find the core with MCS, then decompose the series on it", {tag: ["@guide", "@help:datagrok/solutions/domains/chem"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens sar-small dataset", () => openDataset(page, ds("sar-small")));
    await session.step(17, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
    await session.step(18, "Then \"R-Groups Analysis\" dialog should be visible", () => shouldBe(page, el("\"R-Groups Analysis\" dialog"), "visible"));
    await session.step(19, "When user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")));
    await session.step(20, "Then \"R-Groups Analysis\" dialog should have finished updating", () => finishedUpdating(page, el("\"R-Groups Analysis\" dialog")));
    await session.step(21, "When user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
    await session.step(22, "Then a new column \"R1\" should have been added", () => newColumnNamed(page, "R1"));
    await session.step(23, "And a new column \"R2\" should have been added", () => newColumnNamed(page, "R2"));
    await session.step(24, "And trellis plot viewer should be visible", () => shouldBe(page, el("trellis plot viewer"), "visible"));
  });
});
