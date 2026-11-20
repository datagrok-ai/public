/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/quality/F-009.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {aiPanelNewChat, askAI, liveModelGate, openDemog, runtimeHealthy, screenshot, trafficRecorded} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {propertyShouldBe, propertyShouldNotBe, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A bar chart of averages \"colored by\" a category draws bars, not an error or a silent no-op", () => {
  const session = feature(test, "features/quality/F-009.feature", import.meta.url);
  test("The bar chart of average height by race draws, and is not stacked by disease", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(13, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(14, "And user is logged in", () => loggedIn(page));
    await session.step(15, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(16, "And user opens the demog demo table", () => openDemog(page));
    await session.step(17, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(20, "When user asks the AI panel \"Which race has the highest average height? Show the answer in a bar chart colored by disease, then tell me the top race and its value.\"", () => askAI(page, "Which race has the highest average height? Show the answer in a bar chart colored by disease, then tell me the top race and its value."));
    await session.step(21, "Then the open tableview should have 1 Bar chart viewer", () => viewerCount(page, 1, "Bar chart"));
    await session.step(22, "And \"valueAggrType\" property of bar chart viewer should be \"avg\"", () => propertyShouldBe(page, "valueAggrType", el("bar chart viewer"), "avg"));
    await session.step(23, "And \"colorColumnName\" property of bar chart viewer should not be \"disease\"", () => propertyShouldNotBe(page, "colorColumnName", el("bar chart viewer"), "disease"));
    await session.step(24, "And \"stackColumnName\" property of bar chart viewer should not be \"disease\"", () => propertyShouldNotBe(page, "stackColumnName", el("bar chart viewer"), "disease"));
    await session.step(25, "And a screenshot \"F-009-bar-chart\" of the AI panel is saved", () => screenshot(page, "F-009-bar-chart"));
  });
});
