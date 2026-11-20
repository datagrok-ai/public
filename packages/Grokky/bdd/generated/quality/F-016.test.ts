/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/quality/F-016.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {aiPanelNewChat, askAI, liveModelGate, runtimeHealthy, screenshot, trafficRecorded} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noViewerError, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Asked to fix uniform bar charts of averages stacked by a category, the assistant removes the stack", () => {
  const session = feature(test, "features/quality/F-016.feature", import.meta.url);
  test("The stacks under avg are removed, and the columns keep no color coding", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(12, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(13, "And user is logged in", () => loggedIn(page));
    await session.step(14, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(15, "And user opens a table \"patients\" with:", () => openTableOf(page, "patients", [["study","site","sex","disease","age","weight"],["S1","A","F","Asthma","34","61"],["S1","A","M","Diabetes","52","88"],["S1","B","F","None","45","70"],["S1","B","M","Asthma","29","79"],["S2","A","F","Diabetes","61","74"],["S2","A","M","None","38","92"],["S2","B","F","Asthma","57","66"],["S2","B","M","Diabetes","44","84"],["S3","C","F","None","23","58"],["S3","C","M","Asthma","66","81"]]), [["study","site","sex","disease","age","weight"],["S1","A","F","Asthma","34","61"],["S1","A","M","Diabetes","52","88"],["S1","B","F","None","45","70"],["S1","B","M","Asthma","29","79"],["S2","A","F","Diabetes","61","74"],["S2","A","M","None","38","92"],["S2","B","F","Asthma","57","66"],["S2","B","M","Diabetes","44","84"],["S3","C","F","None","23","58"],["S3","C","M","Asthma","66","81"]]);
    await session.step(27, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","study"],["valueColumnName","weight"],["valueAggrType","avg"],["stackColumnName","sex"]]), [["splitColumnName","study"],["valueColumnName","weight"],["valueAggrType","avg"],["stackColumnName","sex"]]);
    await session.step(32, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","site"],["valueColumnName","age"],["valueAggrType","avg"],["stackColumnName","disease"]]), [["splitColumnName","site"],["valueColumnName","age"],["valueAggrType","avg"],["stackColumnName","disease"]]);
    await session.step(37, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(40, "When user asks the AI panel \"Every bar in both charts looks the same green to me. Can you check what you see and fix it?\"", () => askAI(page, "Every bar in both charts looks the same green to me. Can you check what you see and fix it?"));
    await session.step(41, "Then \"stackColumnName\" property of first bar chart viewer should be \"\"", () => propertyShouldBe(page, "stackColumnName", el("first bar chart viewer"), ""));
    await session.step(42, "And \"stackColumnName\" property of second bar chart viewer should be \"\"", () => propertyShouldBe(page, "stackColumnName", el("second bar chart viewer"), ""));
    await session.step(43, "And \"sex\" column should have no color coding", () => noColorCoding(page, "sex"));
    await session.step(44, "And \"disease\" column should have no color coding", () => noColorCoding(page, "disease"));
    await session.step(45, "And no viewer of the current view should report an error", () => noViewerError(page));
    await session.step(46, "And a screenshot \"F-016-bar-charts\" of the AI panel is saved", () => screenshot(page, "F-016-bar-charts"));
  });
});
