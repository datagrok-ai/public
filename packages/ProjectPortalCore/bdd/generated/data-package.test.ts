/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/data-package.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {exampleRunSaved, openApp, openTableStarting, scratchProjectAtStage} from '../bindings/steps.js';
import {clickOn, shouldBe, shouldContainText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Data package manifest", () => {
  const session = feature(test, "features/data-package.feature", import.meta.url);
  test("A project's data package lists its criteria and its saved run", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Given a scratch project \"BDD-PKG\" at stage \"D5\" in the \"EX-BIOC\" portfolio", () => scratchProjectAtStage(page, "BDD-PKG", "D5", "EX-BIOC"));
    await session.step(13, "And a run of the example model saved to the scratch project \"BDD-PKG\" in the \"EX-BIOC\" portfolio", () => exampleRunSaved(page, "BDD-PKG", "EX-BIOC"));
    await session.step(14, "When user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(15, "And user clicks on scratch package project card", () => clickOn(page, el("scratch package project card")));
    await session.step(16, "And user clicks on data package button", () => clickOn(page, el("data package button")));
    await session.step(17, "Then data package dialog should be visible", () => shouldBe(page, el("data package dialog"), "visible"));
    await session.step(18, "And data package summary should contain text \"1 run\"", () => shouldContainText(page, el("data package summary"), "1 run"));
    await session.step(19, "And data package lines should contain text \"Fed-batch titer (example) · BDD-PKG\"", () => shouldContainText(page, el("data package lines"), "Fed-batch titer (example) · BDD-PKG"));
    await session.step(20, "And data package lines should contain text \"Production bioreactor\"", () => shouldContainText(page, el("data package lines"), "Production bioreactor"));
    await session.step(21, "And data package lines should contain text \"Process description\"", () => shouldContainText(page, el("data package lines"), "Process description"));
    await session.step(22, "And download data package csv button should be visible", () => shouldBe(page, el("download data package csv button"), "visible"));
    await session.step(23, "When user clicks on open data package as table button", () => clickOn(page, el("open data package as table button")));
    await session.step(24, "Then data package dialog should not be visible", () => shouldNotBe(page, el("data package dialog"), "visible"));
    await session.step(25, "And the open table should have 9 rows and a name starting with \"BDD-PKG data package\"", () => openTableStarting(page, 9, "BDD-PKG data package"));
  });
});
