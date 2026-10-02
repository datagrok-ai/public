/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/profile.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp} from '../bindings/steps.js';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A portfolio exported as its file", () => {
  const session = feature(test, "features/profile.feature", import.meta.url);
  test("The board exports the open portfolio", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(12, "And user clicks on export button", () => clickOn(page, el("export button")));
    await session.step(13, "Then export dialog should be visible", () => shouldBe(page, el("export dialog"), "visible"));
    await session.step(14, "And exported portfolio file should contain text \"EX-SMD-V1\"", () => shouldContainText(page, el("exported portfolio file"), "EX-SMD-V1"));
    await session.step(15, "And exported portfolio file should contain text \"Hit-to-Lead\"", () => shouldContainText(page, el("exported portfolio file"), "Hit-to-Lead"));
  });
});
