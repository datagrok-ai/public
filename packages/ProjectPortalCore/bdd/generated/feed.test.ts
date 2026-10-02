/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/feed.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, scratchProjectAtStage} from '../bindings/steps.js';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Activity feed", () => {
  const session = feature(test, "features/feed.feature", import.meta.url);
  test("A gate decision shows in the project's activity and in My feed", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Given a scratch project \"BDD-FEED\" at stage \"TV\" in the \"EX-SMD\" portfolio", () => scratchProjectAtStage(page, "BDD-FEED", "TV", "EX-SMD"));
    await session.step(13, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(14, "And user clicks on feed project card", () => clickOn(page, el("feed project card")));
    await session.step(15, "And user clicks on follow button", () => clickOn(page, el("follow button")));
    await session.step(16, "And user clicks on assay criterion met button", () => clickOn(page, el("assay criterion met button")));
    await session.step(17, "And user clicks on advance button", () => clickOn(page, el("advance button")));
    await session.step(18, "Then activity heading should contain text \"Activity\"", () => shouldContainText(page, el("activity heading"), "Activity"));
    await session.step(19, "And activity entries should contain text \"Advance: Target Validation →\"", () => shouldContainText(page, el("activity entries"), "Advance: Target Validation →"));
    await session.step(20, "When user selects \"Gate decisions (1)\" in activity kind input", () => selectIn(page, "Gate decisions (1)", el("activity kind input")));
    await session.step(21, "Then activity entries should contain text \"with 1 piece of evidence\"", () => shouldContainText(page, el("activity entries"), "with 1 piece of evidence"));
    await session.step(22, "When user clicks on my feed button", () => clickOn(page, el("my feed button")));
    await session.step(23, "Then my feed dialog should be visible", () => shouldBe(page, el("my feed dialog"), "visible"));
    await session.step(24, "And my feed entries should contain text \"Advance: Target Validation →\"", () => shouldContainText(page, el("my feed entries"), "Advance: Target Validation →"));
    await session.step(25, "When user clicks on first my feed program link", () => clickOn(page, el("first my feed program link")));
    await session.step(26, "Then my feed dialog should not be visible", () => shouldNotBe(page, el("my feed dialog"), "visible"));
    await session.step(27, "And program title should contain text \"BDD-FEED\"", () => shouldContainText(page, el("program title"), "BDD-FEED"));
  });
});
