/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/quality/F-004.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {aiPanelNewChat, askAI, barsDrawnInOrder, liveModelGate, openDemog, runtimeHealthy, screenshot, trafficRecorded} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {addCalculated} from '@datagrok-libraries/bdd/bindings/platform/data';
import {addViewerWith} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Asked to order a bar chart's categories, the assistant changes the order the bars are drawn in", () => {
  const session = feature(test, "features/quality/F-004.feature", import.meta.url);
  test("The age groups of a bar chart are drawn youngest to oldest", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(13, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(14, "And user is logged in", () => loggedIn(page));
    await session.step(15, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(16, "And user opens the demog demo table", () => openDemog(page));
    await session.step(17, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(20, "Given user adds a calculated column \"age group\" with formula 'If(${age} < 30, \"under 30\", If(${age} <= 50, \"30 to 50\", \"over 50\"))'", () => addCalculated(page, "age group", "If(${age} < 30, \"under 30\", If(${age} <= 50, \"30 to 50\", \"over 50\"))"));
    await session.step(21, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","age group"],["stackColumnName","sex"],["valueAggrType","count"]]), [["splitColumnName","age group"],["stackColumnName","sex"],["valueAggrType","count"]]);
    await session.step(25, "When user asks the AI panel \"The order of the bars looks odd, can you put the age groups youngest to oldest?\"", () => askAI(page, "The order of the bars looks odd, can you put the age groups youngest to oldest?"));
    await session.step(26, "Then the bar chart should draw its bars in the order \"under 30, 30 to 50, over 50\"", () => barsDrawnInOrder(page, "under 30, 30 to 50, over 50"));
    await session.step(27, "And a screenshot \"F-004-bar-order\" of the AI panel is saved", () => screenshot(page, "F-004-bar-order"));
  });
});
