/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/modes.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, preferencesCleared} from '../bindings/steps.js';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Simple and advanced mode, quick actions per persona", () => {
  const session = feature(test, "features/modes.feature", import.meta.url);
  test("Switch persona and mode on a discovery project and see the page change", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "Given the user's portal mode and persona are not set", () => preferencesCleared(page));
    await session.step(12, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(13, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(14, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(15, "And user selects \"Chemist\" in persona input", () => selectIn(page, "Chemist", el("persona input")));
    await session.step(16, "Then quick actions heading should contain text \"Chemist\"", () => shouldContainText(page, el("quick actions heading"), "Chemist"));
    await session.step(17, "And quick actions should contain text \"Open SAR table\"", () => shouldContainText(page, el("quick actions"), "Open SAR table"));
    await session.step(18, "And quick actions should contain text \"Needs attention\"", () => shouldContainText(page, el("quick actions"), "Needs attention"));
    await session.step(19, "And saved runs heading should be visible", () => shouldBe(page, el("saved runs heading"), "visible"));
    await session.step(20, "When user selects \"Modeler\" in persona input", () => selectIn(page, "Modeler", el("persona input")));
    await session.step(21, "Then quick actions heading should contain text \"Modeler\"", () => shouldContainText(page, el("quick actions heading"), "Modeler"));
    await session.step(22, "And quick actions should contain text \"Save a run\"", () => shouldContainText(page, el("quick actions"), "Save a run"));
    await session.step(23, "And quick actions should not contain text \"Needs attention\"", () => shouldNotContainText(page, el("quick actions"), "Needs attention"));
    await session.step(24, "When user selects \"Simple\" in mode input", () => selectIn(page, "Simple", el("mode input")));
    await session.step(25, "Then saved runs heading should be absent", () => shouldBe(page, el("saved runs heading"), "absent"));
    await session.step(26, "And more button should contain text \"More\"", () => shouldContainText(page, el("more button"), "More"));
    await session.step(27, "When user clicks on more button", () => clickOn(page, el("more button")));
    await session.step(28, "Then saved runs heading should be visible", () => shouldBe(page, el("saved runs heading"), "visible"));
    await session.step(29, "When user selects \"Advanced\" in mode input", () => selectIn(page, "Advanced", el("mode input")));
    await session.step(30, "Then more button should be absent", () => shouldBe(page, el("more button"), "absent"));
  });
});
