/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/gates.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, scratchProject} from '../bindings/steps.js';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Gate decisions follow the portfolio", () => {
  const session = feature(test, "features/gates.feature", import.meta.url);
  test("A registry-owned gate cannot be advanced here", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "When user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(12, "And user clicks on EXB-201 card", () => clickOn(page, el("EXB-201 card")));
    await session.step(13, "Then stage track should be visible", () => shouldBe(page, el("stage track"), "visible"));
    await session.step(14, "And advance button should be disabled", () => shouldBe(page, el("advance button"), "disabled"));
    await session.step(15, "And gate reasons should contain text \"The registry owns this portfolio's stages\"", () => shouldContainText(page, el("gate reasons"), "The registry owns this portfolio's stages"));
  });
  test("A discovery project enters its model and stops at the first gate with criteria", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(18, "Given a scratch project \"BDD-SCRATCH\" in the \"EX-SMD\" portfolio", () => scratchProject(page, "BDD-SCRATCH", "EX-SMD"));
    await session.step(19, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(20, "And user clicks on scratch project card", () => clickOn(page, el("scratch project card")));
    await session.step(21, "Then advance button should contain text \"Target Validation\"", () => shouldContainText(page, el("advance button"), "Target Validation"));
    await session.step(22, "When user clicks on advance button", () => clickOn(page, el("advance button")));
    await session.step(23, "Then scratch project card in target validation column should be visible", () => shouldBe(page, el("scratch project card in target validation column"), "visible"));
    await session.step(24, "And transition history should contain text \"advance\"", () => shouldContainText(page, el("transition history"), "advance"));
    await session.step(25, "And advance button should be disabled", () => shouldBe(page, el("advance button"), "disabled"));
    await session.step(26, "And gate reasons should contain text \"required criteria of Hit Finding\"", () => shouldContainText(page, el("gate reasons"), "required criteria of Hit Finding"));
    await session.step(27, "When user clicks on assay criterion met button", () => clickOn(page, el("assay criterion met button")));
    await session.step(28, "Then advance button should be enabled", () => shouldBe(page, el("advance button"), "enabled"));
  });
});
