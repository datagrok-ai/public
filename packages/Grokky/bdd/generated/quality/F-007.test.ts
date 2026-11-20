/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/quality/F-007.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {aiPanelNewChat, askAI, liveModelGate, nothingAfterStop, runtimeHealthy, screenshot, sendToAI, trafficRecorded, turnStopped} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, waitSeconds} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A request the user stopped stays stopped", () => {
  const session = feature(test, "features/quality/F-007.feature", import.meta.url);
  test("A punctuation-only prompt after Stop does not redo the stopped job", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(12, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(13, "And user is logged in", () => loggedIn(page));
    await session.step(14, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(15, "And user opens mol1K dataset", () => openDataset(page, ds("mol1K")));
    await session.step(16, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(19, "When user sends \"Add 8 columns one by one named stop_1 to stop_8, each equal to Q plus its number, one step each\" to the AI panel", () => sendToAI(page, "Add 8 columns one by one named stop_1 to stop_8, each equal to Q plus its number, one step each"));
    await session.step(20, "And user waits 6 seconds", () => waitSeconds(page, 6));
    await session.step(21, "And user clicks on AI stop button", () => clickOn(page, el("AI stop button")));
    await session.step(22, "Then the AI turn should have been stopped", () => turnStopped(page));
    await session.step(23, "And the table should not have a column \"stop_8\"", () => hasNoColumn(page, "stop_8"));
    await session.step(24, "When user asks the AI panel \"How many rows does the table have?\"", () => askAI(page, "How many rows does the table have?"));
    await session.step(25, "Then the table should not have a column \"stop_8\"", () => hasNoColumn(page, "stop_8"));
    await session.step(26, "When user asks the AI panel \"... ??? !!!\"", () => askAI(page, "... ??? !!!"));
    await session.step(27, "Then the table should not have a column \"stop_8\"", () => hasNoColumn(page, "stop_8"));
    await session.step(28, "And the stopped AI turn should have sent nothing after it was aborted", () => nothingAfterStop(page));
    await session.step(29, "And a screenshot \"F-007-after-punctuation\" of the AI panel is saved", () => screenshot(page, "F-007-after-punctuation"));
  });
});
