/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/program-roles.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noProgramRole, openApp} from '../bindings/steps.js';
import {clickOn, enterInto, selectIn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Program roles", () => {
  const session = feature(test, "features/program-roles.feature", import.meta.url);
  test("An approver makes a user a contributor of a discovery project, then removes the role", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "Given \"admin\" has no role on \"EXS-KIN1\" in the \"EX-SMD\" portfolio", () => noProgramRole(page, "admin", "EXS-KIN1", "EX-SMD"));
    await session.step(12, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(13, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(14, "Then program roles should contain text \"None\"", () => shouldContainText(page, el("program roles"), "None"));
    await session.step(15, "When user enters \"admin\" into role principal input", () => enterInto(page, "admin", el("role principal input")));
    await session.step(16, "And user selects \"Contributor\" in program role input", () => selectIn(page, "Contributor", el("program role input")));
    await session.step(17, "And user clicks on add role button", () => clickOn(page, el("add role button")));
    await session.step(18, "Then program roles should contain text \"Contributor\"", () => shouldContainText(page, el("program roles"), "Contributor"));
    await session.step(19, "And role result should contain text \"workspace folder shared with Admin (edit)\"", () => shouldContainText(page, el("role result"), "workspace folder shared with Admin (edit)"));
    await session.step(20, "When user clicks on remove admin role link", () => clickOn(page, el("remove admin role link")));
    await session.step(21, "Then program roles should contain text \"None\"", () => shouldContainText(page, el("program roles"), "None"));
    await session.step(22, "And role result should contain text \"no program role\"", () => shouldContainText(page, el("role result"), "no program role"));
  });
});
