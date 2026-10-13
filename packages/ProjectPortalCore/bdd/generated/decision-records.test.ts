/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/decision-records.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, scratchProjectAtStage} from '../bindings/steps.js';
import {clickOn, enterInto, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Decisions with evidence where the registry owns the stage", () => {
  const session = feature(test, "features/decision-records.feature", import.meta.url);
  test("Record a gate decision with its evidence without moving the stage, then find it", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "Given a scratch project \"BDD-REC\" at stage \"D2\" in the \"EX-BIOC\" portfolio", () => scratchProjectAtStage(page, "BDD-REC", "D2", "EX-BIOC"));
    await session.step(12, "When user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(13, "And user clicks on scratch record project card", () => clickOn(page, el("scratch record project card")));
    await session.step(14, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(15, "Then record decision button should contain text \"Record decision for\"", () => shouldContainText(page, el("record decision button"), "Record decision for"));
    await session.step(16, "And record note should contain text \"the stage is set by the registry\"", () => shouldContainText(page, el("record note"), "the stage is set by the registry"));
    await session.step(17, "And advance button should be disabled", () => shouldBe(page, el("advance button"), "disabled"));
    await session.step(18, "When user enters \"Comparability shown in the gate review\" into rationale input", () => enterInto(page, "Comparability shown in the gate review", el("rationale input")));
    await session.step(19, "And user selects \"Hold\" in record outcome input", () => selectIn(page, "Hold", el("record outcome input")));
    await session.step(20, "And user clicks on record decision button", () => clickOn(page, el("record decision button")));
    await session.step(21, "Then transition history should contain text \"Hold at the gate to\"", () => shouldContainText(page, el("transition history"), "Hold at the gate to"));
    await session.step(22, "And transition history should contain text \"decision recorded; stage set by registry\"", () => shouldContainText(page, el("transition history"), "decision recorded; stage set by registry"));
    await session.step(23, "When user clicks on board link", () => clickOn(page, el("board link")));
    await session.step(24, "Then scratch record project card in gate D2 column should be visible", () => shouldBe(page, el("scratch record project card in gate D2 column"), "visible"));
    await session.step(25, "When user clicks on decisions button", () => clickOn(page, el("decisions button")));
    await session.step(26, "And user enters \"comparability\" into decision text input", () => enterInto(page, "comparability", el("decision text input")));
    await session.step(27, "Then decision results should contain text \"BDD-REC\"", () => shouldContainText(page, el("decision results"), "BDD-REC"));
  });
});
