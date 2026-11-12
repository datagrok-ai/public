/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/ai-panel/table-work.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {aiPanelNewChat, askAI, askAIAndAnswer, bmiValues, groupSummary, hiddenInGrid, highestMeanNamed, liveModelGate, numericColumnsNamed, oneOfToolsRan, onlyGrid, questionAsked, runtimeHealthy, sameSession, screenshot, submitToAI, toolRan, trafficRecorded} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {columnComplete, columnCount, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorCodedAs, filterIsExactly, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, propertiesShouldBe, propertyShouldBe, showsFilteredRows, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("The AI panel works on the open table", () => {
  const session = feature(test, "features/ai-panel/table-work.feature", import.meta.url);
  test("A small dashboard is built from one request, then changed by a follow-up", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(14, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(15, "And user is logged in", () => loggedIn(page));
    await session.step(16, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(21, "When user asks the AI panel \"Help me build a small dashboard on this table. First remove the rows where HEIGHT or WEIGHT is missing. Then add a column named BMI with each patient's body-mass index (WEIGHT is in kg, HEIGHT in cm). Then add a scatter plot with HEIGHT on the X axis and WEIGHT on the Y axis, colored by SEX, and a box plot of BMI for each RACE. Finally, filter the table to patients aged 50 or older.\"", () => askAI(page, "Help me build a small dashboard on this table. First remove the rows where HEIGHT or WEIGHT is missing. Then add a column named BMI with each patient's body-mass index (WEIGHT is in kg, HEIGHT in cm). Then add a scatter plot with HEIGHT on the X axis and WEIGHT on the Y axis, colored by SEX, and a box plot of BMI for each RACE. Finally, filter the table to patients aged 50 or older."));
    await session.step(22, "Then the table should have 872 rows", () => rowCount(page, 872));
    await session.step(23, "And \"HEIGHT\" column should have no missing values", () => columnComplete(page, "HEIGHT"));
    await session.step(24, "And every value of \"BMI\" column should be the BMI of \"WEIGHT\" in kg and \"HEIGHT\" in cm", () => bmiValues(page, "BMI", "WEIGHT", "HEIGHT"));
    await session.step(25, "And the open tableview should have 1 Scatter plot viewer", () => viewerCount(page, 1, "Scatter plot"));
    await session.step(26, "And properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]), [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]);
    await session.step(30, "And the open tableview should have 1 Box plot viewer", () => viewerCount(page, 1, "Box plot"));
    await session.step(31, "And properties of box plot viewer should be:", () => propertiesShouldBe(page, el("box plot viewer"), [["valueColumnName","BMI"],["category1ColumnName","RACE"]]), [["valueColumnName","BMI"],["category1ColumnName","RACE"]]);
    await session.step(34, "And the filter should pass exactly the rows where \"AGE\" is between 50 and 200", () => filterIsExactly(page, "AGE", 50, 200));
    await session.step(35, "And scatter plot viewer should show every row that passes the filter of its table", () => showsFilteredRows(page, el("scatter plot viewer")));
    await session.step(36, "When user asks the AI panel \"On the scatter plot you just made, color the points by RACE instead of SEX.\"", () => askAI(page, "On the scatter plot you just made, color the points by RACE instead of SEX."));
    await session.step(37, "Then the open tableview should have 1 Scatter plot viewer", () => viewerCount(page, 1, "Scatter plot"));
    await session.step(38, "And properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"],["colorColumnName","RACE"]]), [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"],["colorColumnName","RACE"]]);
    await session.step(42, "And the open tableview should have 1 Box plot viewer", () => viewerCount(page, 1, "Box plot"));
    await session.step(43, "And every turn of the chat should have used the same AI session", () => sameSession(page));
    await session.step(44, "And a screenshot \"table-dashboard\" of the AI panel is saved", () => screenshot(page, "table-dashboard"));
  });
  test("A question about the table is answered with the value computed from it", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(14, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(15, "And user is logged in", () => loggedIn(page));
    await session.step(16, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(47, "When user asks the AI panel \"Look at the women in this table: which race is heaviest on average? Reply with just the race and its mean weight, rounded to one decimal.\"", () => askAI(page, "Look at the women in this table: which race is heaviest on average? Reply with just the race and its mean weight, rounded to one decimal."));
    await session.step(48, "Then the AI answer should name the \"RACE\" with the highest mean \"WEIGHT\" where \"SEX\" is \"F\"", () => highestMeanNamed(page, "RACE", "WEIGHT", "SEX", "F"));
    await session.step(49, "And the AI should have used one of the tools \"datagrok_exec, terminal\"", () => oneOfToolsRan(page, "datagrok_exec, terminal"));
    await session.step(50, "And a screenshot \"table-question\" of the AI panel is saved", () => screenshot(page, "table-question"));
  });
  test("The table view's own commands color-code a column and hide another", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(14, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(15, "And user is logged in", () => loggedIn(page));
    await session.step(16, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(53, "When user asks the AI panel \"Using this table view's own commands, apply linear color coding to the AGE column and hide the CONTROL column from the grid.\"", () => askAI(page, "Using this table view's own commands, apply linear color coding to the AGE column and hide the CONTROL column from the grid."));
    await session.step(54, "Then \"AGE\" column should be color-coded linearly", () => colorCodedAs(page, "AGE", "linearly"));
    await session.step(55, "And the \"CONTROL\" column should be hidden in the grid", () => hiddenInGrid(page, "CONTROL"));
    await session.step(56, "And the table should have a column \"CONTROL\"", () => hasColumn(page, "CONTROL"));
    await session.step(57, "And the AI should have used the \"list_view_functions\" tool", () => toolRan(page, "list_view_functions"));
    await session.step(58, "And the AI should have used the \"call_view_function\" tool", () => toolRan(page, "call_view_function"));
    await session.step(59, "And a screenshot \"table-view-functions\" of the AI panel is saved", () => screenshot(page, "table-view-functions"));
  });
  test("The AI asks which chart to add and adds the one picked", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(14, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(15, "And user is logged in", () => loggedIn(page));
    await session.step(16, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(62, "When user asks the AI panel \"I'd like a chart of AGE in this table, but I'm not sure which kind. Before doing anything, ask me which chart type I want and offer me a few options.\" and picks the option matching \"histogram\" in its question", () => askAIAndAnswer(page, "I'd like a chart of AGE in this table, but I'm not sure which kind. Before doing anything, ask me which chart type I want and offer me a few options.", "histogram"));
    await session.step(63, "Then the AI should have asked a question offering 2 or more options", () => questionAsked(page, 2));
    await session.step(64, "And the AI should have used the \"AskUserQuestion\" tool", () => toolRan(page, "AskUserQuestion"));
    await session.step(65, "And the open tableview should have 1 Histogram viewer", () => viewerCount(page, 1, "Histogram"));
    await session.step(66, "And \"valueColumnName\" property of histogram viewer should be \"AGE\"", () => propertyShouldBe(page, "valueColumnName", el("histogram viewer"), "AGE"));
    await session.step(67, "And a screenshot \"table-ask-user\" of the AI panel is saved", () => screenshot(page, "table-ask-user"));
  });
  test("A request for data the table does not have changes nothing", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(14, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(15, "And user is logged in", () => loggedIn(page));
    await session.step(16, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(70, "When user asks the AI panel \"I'd like to see how cholesterol changes with age in these patients. If this table has no cholesterol data, don't make anything up: tell me which numeric columns it has instead.\"", () => askAI(page, "I'd like to see how cholesterol changes with age in these patients. If this table has no cholesterol data, don't make anything up: tell me which numeric columns it has instead."));
    await session.step(71, "Then the AI answer should name every numeric column of the table", () => numericColumnsNamed(page));
    await session.step(72, "And the current table view should show nothing but its grid", () => onlyGrid(page));
    await session.step(73, "And the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(74, "And the table should have 11 columns", () => columnCount(page, 11));
    await session.step(75, "And a screenshot \"table-missing-column\" of the AI panel is saved", () => screenshot(page, "table-missing-column"));
  });
  test("A per-group summary goes into a new table view", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(14, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(15, "And user is logged in", () => loggedIn(page));
    await session.step(16, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await session.step(78, "When user asks the AI panel \"I need a per-RACE summary of this study: how many subjects are in each group, and their average AGE and WEIGHT. Put that summary into a new table and open it in its own view.\"", () => askAI(page, "I need a per-RACE summary of this study: how many subjects are in each group, and their average AGE and WEIGHT. Put that summary into a new table and open it in its own view."));
    await session.step(79, "Then the current table should hold the count and the mean \"AGE, WEIGHT\" of every \"RACE\" group of table \"demog-1000\"", () => groupSummary(page, "AGE, WEIGHT", "RACE", "demog-1000"));
    await session.step(80, "And a screenshot \"table-summary\" of the AI panel is saved", () => screenshot(page, "table-summary"));
  });
  test("A question naming an average and a column is answered, not turned into a chart", {tag: ["@serial", "@known-failure:87"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(14, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(15, "And user is logged in", () => loggedIn(page));
    await session.step(16, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await knownFailure(async () => {
      await session.step(85, "When user submits \"Among the female patients in this table, which RACE has the highest average WEIGHT, and what is that average? Reply with just the race and the value rounded to one decimal.\" in the AI panel", () => submitToAI(page, "Among the female patients in this table, which RACE has the highest average WEIGHT, and what is that average? Reply with just the race and the value rounded to one decimal."));
      await session.step(86, "Then a screenshot \"table-question-intercepted\" of the AI panel is saved", () => screenshot(page, "table-question-intercepted"));
      await session.step(87, "And the AI answer should name the \"RACE\" with the highest mean \"WEIGHT\" where \"SEX\" is \"F\"", () => highestMeanNamed(page, "RACE", "WEIGHT", "SEX", "F"));
    }, 87);
  });
  test("A request to recolor the scatter plot changes it rather than adding another", {tag: ["@serial", "@known-failure:98"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given live AI model calls are allowed", () => liveModelGate(page));
    await session.step(14, "And the AI runtime traffic is recorded", () => trafficRecorded(page));
    await session.step(15, "And user is logged in", () => loggedIn(page));
    await session.step(16, "And the AI runtime container is healthy", () => runtimeHealthy(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And the AI panel is open on a new chat", () => aiPanelNewChat(page));
    await knownFailure(async () => {
      await session.step(92, "Given user adds a Scatter plot viewer with:", () => addViewerWith(page, "Scatter plot", [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]), [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]);
      await session.step(96, "When user submits \"Now recolor that scatter plot by RACE.\" in the AI panel", () => submitToAI(page, "Now recolor that scatter plot by RACE."));
      await session.step(97, "Then a screenshot \"table-recolor-intercepted\" of the AI panel is saved", () => screenshot(page, "table-recolor-intercepted"));
      await session.step(98, "And the open tableview should have 1 Scatter plot viewer", () => viewerCount(page, 1, "Scatter plot"));
      await session.step(99, "And \"colorColumnName\" property of scatter plot viewer should be \"RACE\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "RACE"));
    }, 98);
  });
});
