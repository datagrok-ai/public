/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/follow.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {notFollowed, openApp} from '../bindings/steps.js';
import {clickOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Follow a program and hear about new data", () => {
  const session = feature(test, "features/follow.feature", import.meta.url);
  test("Follow a discovery project and one of its series, then stop", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Given the user does not follow \"EXS-KIN1\" in the \"EX-SMD\" portfolio", () => notFollowed(page, "EXS-KIN1", "EX-SMD"));
    await session.step(13, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(14, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(15, "Then follow button should contain text \"Follow\"", () => shouldContainText(page, el("follow button"), "Follow"));
    await session.step(16, "When user clicks on follow button", () => clickOn(page, el("follow button")));
    await session.step(17, "Then follow button should contain text \"Unfollow\"", () => shouldContainText(page, el("follow button"), "Unfollow"));
    await session.step(18, "And what's new should contain text \"Nothing new since your last visit\"", () => shouldContainText(page, el("what's new"), "Nothing new since your last visit"));
    await session.step(19, "When user clicks on KIN-A follow link", () => clickOn(page, el("KIN-A follow link")));
    await session.step(20, "Then KIN-A follow link should contain text \"Unfollow\"", () => shouldContainText(page, el("KIN-A follow link"), "Unfollow"));
    await session.step(21, "When user clicks on follow button", () => clickOn(page, el("follow button")));
    await session.step(22, "Then follow button should contain text \"Follow\"", () => shouldContainText(page, el("follow button"), "Follow"));
    await session.step(23, "And KIN-A follow link should contain text \"Unfollow\"", () => shouldContainText(page, el("KIN-A follow link"), "Unfollow"));
  });
});
