/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/observations.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, scratchProjectAtStage} from '../bindings/steps.js';
import {clickOn, enterInto, shouldBe, shouldContainText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Observations", () => {
  const session = feature(test, "features/observations.feature", import.meta.url);
  test("An observation on a project is recorded, shows in the activity, and is retracted", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Given a scratch project \"BDD-OBS\" at stage \"TV\" in the \"EX-SMD\" portfolio", () => scratchProjectAtStage(page, "BDD-OBS", "TV", "EX-SMD"));
    await session.step(13, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(14, "And user clicks on scratch observation project card", () => clickOn(page, el("scratch observation project card")));
    await session.step(15, "And user clicks on add observation button", () => clickOn(page, el("add observation button")));
    await session.step(16, "Then observation dialog should be visible", () => shouldBe(page, el("observation dialog"), "visible"));
    await session.step(17, "When user enters \"Selectivity is the risk\" into observation title input", () => enterInto(page, "Selectivity is the risk", el("observation title input")));
    await session.step(18, "And user enters \"Two off-target kinases within 10x in the panel\" into observation finding input", () => enterInto(page, "Two off-target kinases within 10x in the panel", el("observation finding input")));
    await session.step(19, "And user clicks on record observation button", () => clickOn(page, el("record observation button")));
    await session.step(20, "Then observation dialog should not be visible", () => shouldNotBe(page, el("observation dialog"), "visible"));
    await session.step(21, "And observations heading should contain text \"Observations (1)\"", () => shouldContainText(page, el("observations heading"), "Observations (1)"));
    await session.step(22, "And observations should contain text \"the program · Selectivity is the risk\"", () => shouldContainText(page, el("observations"), "the program · Selectivity is the risk"));
    await session.step(23, "And activity entries should contain text \"Observation: Selectivity is the risk\"", () => shouldContainText(page, el("activity entries"), "Observation: Selectivity is the risk"));
    await session.step(24, "When user clicks on first retract observation link", () => clickOn(page, el("first retract observation link")));
    await session.step(25, "And user enters \"The panel was rerun\" into retraction input", () => enterInto(page, "The panel was rerun", el("retraction input")));
    await session.step(26, "And user clicks on retract button", () => clickOn(page, el("retract button")));
    await session.step(27, "Then observations heading should contain text \"Observations (0)\"", () => shouldContainText(page, el("observations heading"), "Observations (0)"));
    await session.step(28, "And observations should contain text \"Retracted\"", () => shouldContainText(page, el("observations"), "Retracted"));
    await session.step(29, "And observations should contain text \"The panel was rerun\"", () => shouldContainText(page, el("observations"), "The panel was rerun"));
  });
});
