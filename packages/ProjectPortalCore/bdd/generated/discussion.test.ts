/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/discussion.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {discussionRemoved, openApp, scratchProjectAtStage, scratchWorkspace} from '../bindings/steps.js';
import {clickOn, enterInto, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Discussion and Ask SME", () => {
  const session = feature(test, "features/discussion.feature", import.meta.url);
  test("A post shows in the program's discussion, and Ask SME shows the context it sends", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "Given a scratch project \"BDD-DISC\" at stage \"H2L\" in the \"EX-SMD\" portfolio", () => scratchProjectAtStage(page, "BDD-DISC", "H2L", "EX-SMD"));
    await session.step(12, "And the workspace of the scratch project \"BDD-DISC\" in the \"EX-SMD\" portfolio is built", () => scratchWorkspace(page, "BDD-DISC", "EX-SMD"));
    await session.step(13, "And the discussion of the scratch project \"BDD-DISC\" in the \"EX-SMD\" portfolio is removed at the end", () => discussionRemoved(page, "BDD-DISC", "EX-SMD"));
    await session.step(14, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(15, "And user clicks on scratch discussion project card", () => clickOn(page, el("scratch discussion project card")));
    await session.step(16, "And user enters \"Is the K1 panel enough for the gate?\" into discussion input", () => enterInto(page, "Is the K1 panel enough for the gate?", el("discussion input")));
    await session.step(17, "And user clicks on discussion post button", () => clickOn(page, el("discussion post button")));
    await session.step(18, "Then discussion heading should contain text \"Discussion (1)\"", () => shouldContainText(page, el("discussion heading"), "Discussion (1)"));
    await session.step(19, "And discussion threads should contain text \"Is the K1 panel enough for the gate?\"", () => shouldContainText(page, el("discussion threads"), "Is the K1 panel enough for the gate?"));
    await session.step(20, "When user clicks on ask SME button", () => clickOn(page, el("ask SME button")));
    await session.step(21, "Then ask SME dialog should be visible", () => shouldBe(page, el("ask SME dialog"), "visible"));
    await session.step(22, "And SME context should contain text \"Program: BDD-DISC\"", () => shouldContainText(page, el("SME context"), "Program: BDD-DISC"));
    await session.step(23, "And SME context should contain text \"Next gate: Lead Optimization\"", () => shouldContainText(page, el("SME context"), "Next gate: Lead Optimization"));
    await session.step(24, "And SME note should contain text \"Example: an EX-SMD approver\"", () => shouldContainText(page, el("SME note"), "Example: an EX-SMD approver"));
  });
});
