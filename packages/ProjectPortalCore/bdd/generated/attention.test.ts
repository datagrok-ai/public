/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/attention.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, openTableCurrent} from '../bindings/steps.js';
import {clickOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Compounds needing attention", () => {
  const session = feature(test, "features/attention.feature", import.meta.url);
  test("The compounds of a discovery project that need attention, and one opened", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(13, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(14, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(15, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(16, "Then needs attention heading should contain text \"Needs attention (4)\"", () => shouldContainText(page, el("needs attention heading"), "Needs attention (4)"));
    await session.step(17, "And hERG not measured count should contain text \"4 Potent, hERG not measured\"", () => shouldContainText(page, el("hERG not measured count"), "4 Potent, hERG not measured"));
    await session.step(18, "And lead without solubility count should contain text \"1 Lead without solubility\"", () => shouldContainText(page, el("lead without solubility count"), "1 Lead without solubility"));
    await session.step(19, "And EXC-1006 attention item should contain text \"no solubility\"", () => shouldContainText(page, el("EXC-1006 attention item"), "no solubility"));
    await session.step(20, "When user clicks on EXC-1006 attention link", () => clickOn(page, el("EXC-1006 attention link")));
    await session.step(21, "Then the open table's current row should have \"EXC-1006\" in \"compound\"", () => openTableCurrent(page, "EXC-1006", "compound"));
  });
});
