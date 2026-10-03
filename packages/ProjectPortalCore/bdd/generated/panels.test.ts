/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/panels.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {exampleDesignSets, openApp, preferencesCleared} from '../bindings/steps.js';
import {clickOn, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Program-page panels per persona", () => {
  const session = feature(test, "features/panels.feature", import.meta.url);
  test("The biologist's page on EXS-KIN1 leads with DMTA; the gate is behind More", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Given the user's portal mode and persona are not set", () => preferencesCleared(page));
    await session.step(13, "And the example design sets are loaded", () => exampleDesignSets(page));
    await session.step(14, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(15, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(16, "And user selects \"Biologist\" in persona input", () => selectIn(page, "Biologist", el("persona input")));
    await session.step(17, "Then panels note should contain text \"Panels for Biologist, from EX-SMD\"", () => shouldContainText(page, el("panels note"), "Panels for Biologist, from EX-SMD"));
    await session.step(18, "And dmta heading should be visible", () => shouldBe(page, el("dmta heading"), "visible"));
    await session.step(19, "And observations heading should be visible", () => shouldBe(page, el("observations heading"), "visible"));
    await session.step(20, "And gate heading should be absent", () => shouldBe(page, el("gate heading"), "absent"));
    await session.step(21, "When user clicks on more button", () => clickOn(page, el("more button")));
    await session.step(22, "Then gate heading should be visible", () => shouldBe(page, el("gate heading"), "visible"));
    await session.step(23, "When user selects \"Chemist\" in persona input", () => selectIn(page, "Chemist", el("persona input")));
    await session.step(24, "Then panels note should be absent", () => shouldBe(page, el("panels note"), "absent"));
    await session.step(25, "And saved runs heading should be visible", () => shouldBe(page, el("saved runs heading"), "visible"));
  });
});
