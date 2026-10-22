/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/ai-panel/analysis-work.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {aiPanelNewChat, askAI, groupSummary, liveModelGate, runtimeHealthy, sameSession, screenshot, trafficRecorded} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {everyValueBetween, everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pivotedAggregationMatches} from '@datagrok-libraries/bdd/bindings/tiers/viewers/pivot-table';
import {propertiesShouldBe, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The AI panel carries out longer analyses on the open table", () => {
  const session = feature(test, "features/ai-panel/analysis-work.feature", import.meta.url);
  test("A histogram is refined over three turns, every change landing on the same viewer", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(17, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(18, "And user is logged in", () => loggedIn(page));
    await session.step(19, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(24, "When user asks the AI panel \"Show the patients' ages as a histogram.\"", () => askAI(page, "Show the patients' ages as a histogram."));
    await session.step(25, "Then the open tableview should have 1 Histogram viewer", () => viewerCount(page, 1, "Histogram"));
    await session.step(26, "And properties of histogram viewer should be:", () => propertiesShouldBe(page, el("histogram viewer"), [["valueColumnName","AGE"]]), [["valueColumnName","AGE"]]);
    await session.step(28, "When user asks the AI panel \"Break that histogram down by SEX, so I can compare men and women.\"", () => askAI(page, "Break that histogram down by SEX, so I can compare men and women."));
    await session.step(29, "Then the open tableview should have 1 Histogram viewer", () => viewerCount(page, 1, "Histogram"));
    await session.step(30, "And properties of histogram viewer should be:", () => propertiesShouldBe(page, el("histogram viewer"), [["valueColumnName","AGE"],["splitColumnName","SEX"]]), [["valueColumnName","AGE"],["splitColumnName","SEX"]]);
    await session.step(33, "When user asks the AI panel \"I'd like finer detail on that chart: twenty bins rather than the default.\"", () => askAI(page, "I'd like finer detail on that chart: twenty bins rather than the default."));
    await session.step(34, "Then the open tableview should have 1 Histogram viewer", () => viewerCount(page, 1, "Histogram"));
    await session.step(35, "And properties of histogram viewer should be:", () => propertiesShouldBe(page, el("histogram viewer"), [["valueColumnName","AGE"],["splitColumnName","SEX"],["bins","20"]]), [["valueColumnName","AGE"],["splitColumnName","SEX"],["bins","20"]]);
    await session.step(39, "And every turn of the chat should have used the same AI session", () => sameSession(page));
    await session.step(40, "And a screenshot \"analysis-histogram\" of the AI panel is saved", () => screenshot(page, "analysis-histogram"));
  });
  test("Rows are filtered, a subgroup selected, and the selection extracted into a table of its own", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(17, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(18, "And user is logged in", () => loggedIn(page));
    await session.step(19, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(43, "When user asks the AI panel \"In this table, show only the patients aged 30 to 50, hiding the others rather than deleting them. Among the patients shown, select the Caucasian women. Then put the selected rows into a new table and open it.\"", () => askAI(page, "In this table, show only the patients aged 30 to 50, hiding the others rather than deleting them. Among the patients shown, select the Caucasian women. Then put the selected rows into a new table and open it."));
    await session.step(44, "Then the table should have 225 rows", () => rowCount(page, 225));
    await session.step(45, "And every value of \"AGE\" column should lie between 30 and 50", () => everyValueBetween(page, "AGE", 30, 50));
    await session.step(46, "And every value of \"SEX\" column should match \"F\"", () => everyValueMatches(page, "SEX", "F"));
    await session.step(47, "And every value of \"RACE\" column should match \"Caucasian\"", () => everyValueMatches(page, "RACE", "Caucasian"));
    await session.step(48, "And a screenshot \"analysis-extract\" of the AI panel is saved", () => screenshot(page, "analysis-extract"));
  });
  test("A pivot table gives the mean weight of every race, one column per sex", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(17, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(18, "And user is logged in", () => loggedIn(page));
    await session.step(19, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(51, "When user asks the AI panel \"Lay out the average WEIGHT as a cross-tab of RACE against SEX in a pivot table.\"", () => askAI(page, "Lay out the average WEIGHT as a cross-tab of RACE against SEX in a pivot table."));
    await session.step(52, "Then the open tableview should have 1 Pivot table viewer", () => viewerCount(page, 1, "Pivot table"));
    await session.step(53, "And the aggregated values of pivot table viewer should match \"avg(WEIGHT)\" grouped by \"RACE\" pivoted on \"SEX\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "avg(WEIGHT)", "RACE", "SEX"));
    await session.step(54, "And a screenshot \"analysis-pivot\" of the AI panel is saved", () => screenshot(page, "analysis-pivot"));
  });
  test("A per-race summary is charted in its own view", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(17, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(18, "And user is logged in", () => loggedIn(page));
    await session.step(19, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(57, "When user asks the AI panel \"Create a per-RACE summary table (patient count and average WEIGHT) as a new table in its own view, then chart that summary's average WEIGHT per RACE as bars in the new view.\"", () => askAI(page, "Create a per-RACE summary table (patient count and average WEIGHT) as a new table in its own view, then chart that summary's average WEIGHT per RACE as bars in the new view."));
    await session.step(58, "Then the current table should hold the count and the mean \"WEIGHT\" of every \"RACE\" group of table \"demog-1000\"", () => groupSummary(page, "WEIGHT", "RACE", "demog-1000"));
    await session.step(59, "And the open tableview should have 1 Bar chart viewer", () => viewerCount(page, 1, "Bar chart"));
    await session.step(60, "And properties of bar chart viewer should be:", () => propertiesShouldBe(page, el("bar chart viewer"), [["splitColumnName","RACE"]]), [["splitColumnName","RACE"]]);
    await session.step(62, "And a screenshot \"analysis-summary-chart\" of the AI panel is saved", () => screenshot(page, "analysis-summary-chart"));
  });
});
