/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/quality/F-012.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {aiPanelNewChat, askAI, categoriesDrawnInOrder, columnCategoryOrder, liveModelGate, runtimeHealthy, screenshot, trafficRecorded} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Ordering a new viewer of a column leaves the other viewers of that column in their order", () => {
  const session = feature(test, "features/quality/F-012.feature", import.meta.url);
  test("A box plot ordered light first keeps the bar chart's heavy-first order", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(14, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(15, "And user is logged in", () => loggedIn(page));
    await session.step(16, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(17, "And user opens a table \"patients\" with:", () => openTableOf(page, "patients", [["weight class","age"],["light","25"],["light","31"],["light","44"],["light","52"],["heavy","38"],["heavy","61"],["heavy","47"],["medium","29"],["medium","55"]]), [["weight class","age"],["light","25"],["light","31"],["light","44"],["light","52"],["heavy","38"],["heavy","61"],["heavy","47"],["medium","29"],["medium","55"]]);
    await session.step(28, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(31, "Given the categories of \"weight class\" column are ordered \"heavy, medium, light\"", () => columnCategoryOrder(page, "weight class", "heavy, medium, light"));
    await session.step(32, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","weight class"],["valueColumnName","age"],["valueAggrType","count"],["barSortType","custom"],["barSortOrder","asc"]]), [["splitColumnName","weight class"],["valueColumnName","age"],["valueAggrType","count"],["barSortType","custom"],["barSortOrder","asc"]]);
    await session.step(38, "And the \"Bar chart\" viewer should draw its categories in the order \"heavy, medium, light\"", () => categoriesDrawnInOrder(page, "Bar chart", "heavy, medium, light"));
    await session.step(39, "When user asks the AI panel \"Now also a box plot of age by weight class, light first, then medium, then heavy.\"", () => askAI(page, "Now also a box plot of age by weight class, light first, then medium, then heavy."));
    await session.step(40, "Then the \"Box plot\" viewer should draw its categories in the order \"light, medium, heavy\"", () => categoriesDrawnInOrder(page, "Box plot", "light, medium, heavy"));
    await session.step(41, "And the \"Bar chart\" viewer should draw its categories in the order \"heavy, medium, light\"", () => categoriesDrawnInOrder(page, "Bar chart", "heavy, medium, light"));
    await session.step(42, "And a screenshot \"F-012-box-and-bar-order\" of the AI panel is saved", () => screenshot(page, "F-012-box-and-bar-order"));
  });
});
