/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/ai-panel/session.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {agentsKilled, aiPanelNewChat, askAI, earlierTurnsShown, historyKept, hostNameAnswered, liveModelGate, medianDivided, medianGiven, notRunningInRuntime, nothingAfterStop, oneOfToolsRan, promptMode, reopenChat, runningInRuntime, runsHermes, runtimeHealthy, sameSession, screenshot, sendToAI, systemMessage, toolRan, trafficRecorded, turnStopped} from '../../bindings/steps.js';
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The AI panel's chat survives a stop, a bash command, a reload and a lost agent process", () => {
  const session = feature(test, "features/ai-panel/session.feature", import.meta.url);
  test("A long job stopped mid-run leaves nothing running, and the chat goes on", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(15, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(16, "And user is logged in", () => loggedIn(page));
    await session.step(17, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(22, "When user sends \"In your terminal, run `sleep 613 && echo finished` (it stands in for a long data job) and tell me what it printed.\" to the AI panel", () => sendToAI(page, "In your terminal, run `sleep 613 && echo finished` (it stands in for a long data job) and tell me what it printed."));
    await session.step(23, "Then the AI runtime container should be running \"sleep 613\"", () => runningInRuntime(page, "sleep 613"));
    await session.step(24, "When user clicks on AI stop button", () => clickOn(page, el("AI stop button")));
    await session.step(25, "Then the AI turn should have been stopped", () => turnStopped(page));
    await session.step(26, "And the AI panel should show the message \"Processing aborted by user\"", () => systemMessage(page, "Processing aborted by user"));
    await session.step(27, "And the AI runtime container should not be running \"sleep 613\"", () => notRunningInRuntime(page, "sleep 613"));
    await session.step(28, "When user asks the AI panel \"What is the median WEIGHT of the patients in this table, to one decimal? Reply with just the number.\"", () => askAI(page, "What is the median WEIGHT of the patients in this table, to one decimal? Reply with just the number."));
    await session.step(29, "Then the AI answer should be the median of \"WEIGHT\"", () => medianGiven(page, "WEIGHT"));
    await session.step(30, "And the AI should have used one of the tools \"datagrok_exec, terminal\"", () => oneOfToolsRan(page, "datagrok_exec, terminal"));
    await session.step(31, "And the stopped AI turn should have sent nothing after it was aborted", () => nothingAfterStop(page));
    await session.step(32, "And every turn of the chat should have used the same AI session", () => sameSession(page));
    await session.step(33, "And the AI runtime should have kept the chat's history itself", () => historyKept(page));
    await session.step(34, "And a screenshot \"session-stop\" of the AI panel is saved", () => screenshot(page, "session-stop"));
  });
  test("A !command shows its raw output, and the next prompt is back in full mode", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(15, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(16, "And user is logged in", () => loggedIn(page));
    await session.step(17, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(37, "When user asks the AI panel \"!hostname\"", () => askAI(page, "!hostname"));
    await session.step(38, "Then the AI answer should be exactly the host name of the AI runtime container", () => hostNameAnswered(page));
    await session.step(39, "And the AI should have used the \"terminal\" tool", () => toolRan(page, "terminal"));
    await session.step(40, "And the last prompt should have reached the runtime in bash mode", () => promptMode(page, "bash"));
    await session.step(41, "When user asks the AI panel \"What is the median AGE of the patients in this table, to one decimal? Reply with just the number.\"", () => askAI(page, "What is the median AGE of the patients in this table, to one decimal? Reply with just the number."));
    await session.step(42, "Then the AI answer should be the median of \"AGE\"", () => medianGiven(page, "AGE"));
    await session.step(43, "And the last prompt should have reached the runtime in full mode", () => promptMode(page, "full"));
    await session.step(44, "And every turn of the chat should have used the same AI session", () => sameSession(page));
    await session.step(45, "And a screenshot \"session-bash\" of the AI panel is saved", () => screenshot(page, "session-bash"));
  });
  test("The chat goes on after the browser reloads", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(15, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(16, "And user is logged in", () => loggedIn(page));
    await session.step(17, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(48, "When user asks the AI panel \"What is the median HEIGHT of the patients in this table, in cm, to one decimal? Reply with just the number.\"", () => askAI(page, "What is the median HEIGHT of the patients in this table, in cm, to one decimal? Reply with just the number."));
    await session.step(49, "Then the AI answer should be the median of \"HEIGHT\"", () => medianGiven(page, "HEIGHT"));
    await session.step(50, "When user reloads the page", () => reloadPage(page));
    await session.step(51, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(52, "And user reopens the AI chat from the chat history", () => reopenChat(page));
    await session.step(53, "Then the AI panel should show the earlier turns of the chat", () => earlierTurnsShown(page));
    await session.step(54, "When user asks the AI panel \"What is that in inches, to one decimal? Reply with just the number.\"", () => askAI(page, "What is that in inches, to one decimal? Reply with just the number."));
    await session.step(55, "Then the AI answer should be the median of \"HEIGHT\" divided by 2.54", () => medianDivided(page, "HEIGHT", 2.54));
    await session.step(56, "And every turn of the chat should have used the same AI session", () => sameSession(page));
    await session.step(57, "And the AI runtime should have kept the chat's history itself", () => historyKept(page));
    await session.step(58, "And a screenshot \"session-reload\" of the AI panel is saved", () => screenshot(page, "session-reload"));
  });
  test("The chat goes on after the runtime loses the chat's agent process", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(15, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(16, "And user is logged in", () => loggedIn(page));
    await session.step(17, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(61, "When user asks the AI panel \"What is the median HEIGHT of the patients in this table, in cm, to one decimal? Reply with just the number.\"", () => askAI(page, "What is the median HEIGHT of the patients in this table, in cm, to one decimal? Reply with just the number."));
    await session.step(62, "Then the AI answer should be the median of \"HEIGHT\"", () => medianGiven(page, "HEIGHT"));
    await session.step(63, "And the AI runtime container should run the Hermes agent and not the Claude Agent SDK", () => runsHermes(page));
    await session.step(64, "When the AI runtime's agent processes are killed", () => agentsKilled(page));
    await session.step(65, "And user asks the AI panel \"What is that in inches, to one decimal? Reply with just the number.\"", () => askAI(page, "What is that in inches, to one decimal? Reply with just the number."));
    await session.step(66, "Then the AI answer should be the median of \"HEIGHT\" divided by 2.54", () => medianDivided(page, "HEIGHT", 2.54));
    await session.step(67, "And every turn of the chat should have used the same AI session", () => sameSession(page));
    await session.step(68, "And the AI runtime should have kept the chat's history itself", () => historyKept(page));
    await session.step(69, "And a screenshot \"session-agent-killed\" of the AI panel is saved", () => screenshot(page, "session-agent-killed"));
  });
});
