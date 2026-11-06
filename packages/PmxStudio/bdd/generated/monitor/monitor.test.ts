/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/monitor/monitor.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {finishedRunning} from '../../bindings/pmx.js';
import {openMonitor} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Monitor a running study", () => {
  const session = feature(test, "features/monitor/monitor.feature", import.meta.url);
  test("Follow the example study and reveal the groups", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(6, "Given user is logged in", () => loggedIn(page));
    await session.step(7, "And user opens PMX Studio on Monitor", () => openMonitor(page));
    await session.step(11, "When user clicks on \"Follow an example study\" button", () => clickOn(page, el("\"Follow an example study\" button")));
    await session.step(12, "Then the study should have finished running", () => finishedRunning(page));
    await session.step(13, "And study error should be hidden", () => shouldBe(page, el("study error"), "hidden"));
    await session.step(14, "And \"Needs attention\" section should be visible", () => shouldBe(page, el("\"Needs attention\" section"), "visible"));
    await session.step(15, "And \"Measurement sheet\" section should be visible", () => shouldBe(page, el("\"Measurement sheet\" section"), "visible"));
    await session.step(16, "When user clicks on \"Decisions\" tab", () => clickOn(page, el("\"Decisions\" tab")));
    await session.step(17, "And user clicks on \"Reveal the groups\" button", () => clickOn(page, el("\"Reveal the groups\" button")));
    await session.step(18, "Then \"Log\" section should contain text \"Groups revealed\"", () => shouldContainText(page, el("\"Log\" section"), "Groups revealed"));
  });
});
