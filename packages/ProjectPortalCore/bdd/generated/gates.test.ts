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
  test("Data conditions are evaluated from the program's data", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(31, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(32, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(33, "And user clicks on H2L step", () => clickOn(page, el("H2L step")));
    await session.step(34, "Then gate heading should contain text \"Hit-to-Lead\"", () => shouldContainText(page, el("gate heading"), "Hit-to-Lead"));
    await session.step(35, "And readiness score should contain text \"2 of 2 required\"", () => shouldContainText(page, el("readiness score"), "2 of 2 required"));
    await session.step(36, "And potency criterion should contain text \"from data\"", () => shouldContainText(page, el("potency criterion"), "from data"));
    await session.step(37, "When user clicks on LO step", () => clickOn(page, el("LO step")));
    await session.step(38, "Then gate heading should contain text \"Lead Optimization\"", () => shouldContainText(page, el("gate heading"), "Lead Optimization"));
    await session.step(39, "And readiness score should contain text \"0 of 1 required\"", () => shouldContainText(page, el("readiness score"), "0 of 1 required"));
    await session.step(40, "And ADME panel criterion detail should contain text \"EXC-1006 lacks\"", () => shouldContainText(page, el("ADME panel criterion detail"), "EXC-1006 lacks"));
  });
});
