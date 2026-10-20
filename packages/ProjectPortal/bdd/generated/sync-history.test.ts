/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sync-history.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, portfolioSynced} from '../bindings/steps.js';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A portfolio's sync history", () => {
  const session = feature(test, "features/sync-history.feature", import.meta.url);
  test("The sync history of a portfolio with a schedule", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Given a sync of the \"EX-SMD\" portfolio has run", () => portfolioSynced(page, "EX-SMD"));
    await session.step(13, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(14, "And user clicks on sync history button", () => clickOn(page, el("sync history button")));
    await session.step(15, "Then sync history dialog should be visible", () => shouldBe(page, el("sync history dialog"), "visible"));
    await session.step(16, "And sync schedule should contain text \"Scheduled: every 6h\"", () => shouldContainText(page, el("sync schedule"), "Scheduled: every 6h"));
    await session.step(17, "And sync owners should contain text \"admin\"", () => shouldContainText(page, el("sync owners"), "admin"));
    await session.step(18, "And latest sync run should contain text \"by hand\"", () => shouldContainText(page, el("latest sync run"), "by hand"));
  });
  test("Data advances stages only as far as the portfolio allows", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(21, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(22, "And user clicks on sync history button", () => clickOn(page, el("sync history button")));
    await session.step(23, "Then sync auto-advance should contain text \"Stages advance by data through Hit-to-Lead\"", () => shouldContainText(page, el("sync auto-advance"), "Stages advance by data through Hit-to-Lead"));
    await session.step(24, "When user clicks on preview advance by data button", () => clickOn(page, el("preview advance by data button")));
    await session.step(25, "Then advance by data preview should contain text \"would advance by data (through H2L)\"", () => shouldContainText(page, el("advance by data preview"), "would advance by data (through H2L)"));
    await session.step(26, "And advance by data preview should contain text \"EXS-KIN1 stays: Lead Optimization is past Hit-to-Lead\"", () => shouldContainText(page, el("advance by data preview"), "EXS-KIN1 stays: Lead Optimization is past Hit-to-Lead"));
  });
});
