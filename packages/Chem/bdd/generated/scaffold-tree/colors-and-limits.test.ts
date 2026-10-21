/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scaffold-tree/colors-and-limits.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.scaffold-tree-add-filter]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {treeBuilt} from '../../bindings/scaffold-tree.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {distinctValues, hasColumn, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCodedCategorically} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openDatasetRowsAs, standRunsService, switchTableView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, clickArea, hoverArea, legendLists, noErrors, painted, propertyShouldBe, readingDoesNotRead, readingIs, readingNotAsRemembered, readingReads, rememberReading, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector, readingIncludes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scaffold Tree colors, blocked generation, large tables and two tables", () => {
  const session = feature(test, "features/scaffold-tree/colors-and-limits.feature", import.meta.url);
  test("Scaffold Tree colors, blocked generation, large tables and two tables", {tag: ["@journey", "@realizes:chem.cp.scaffold-tree-add-filter"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(17, "And user picks \"Chem > Analyze > Scaffold Tree\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Scaffold Tree"));
    await run.scenario("The magic wand is blocked on a table with no molecule column", async () => {
      await session.step(20, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
      await session.step(21, "And user adds Scaffold Tree viewer", () => addViewer(page, "Scaffold Tree"));
      await session.step(22, "Then Scaffold Tree viewer should be visible", () => shouldBe(page, el("Scaffold Tree viewer"), "visible"));
      await session.step(23, "And the \"generate blocked reason\" reading of Scaffold Tree viewer should be \"There is no molecule column in the table\"", () => readingReads(page, "generate blocked reason", el("Scaffold Tree viewer"), "There is no molecule column in the table"));
      await session.step(24, "And the \"message\" reading of Scaffold Tree viewer should include the text \"No molecule column found\"", () => readingIncludes(page, "message", el("Scaffold Tree viewer"), "No molecule column found"));
      await session.step(25, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(26, "Then \"Generate\" icon inside Scaffold Tree viewer should be disabled", () => shouldBe(page, el("\"Generate\" icon inside Scaffold Tree viewer"), "disabled"));
      await session.step(27, "And the \"nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "nodes", el("Scaffold Tree viewer"), 0));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On 1,000 molecules the magic wand is open, and says it generates from a sample of 750", async () => {
      await session.step(31, "Given user opens mol1K dataset", () => openDataset(page, ds("mol1K")));
      await session.step(32, "When user picks \"Chem > Analyze > Scaffold Tree\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Scaffold Tree"));
      await session.step(33, "Then Scaffold Tree viewer should be visible", () => shouldBe(page, el("Scaffold Tree viewer"), "visible"));
      await session.step(34, "And the \"generate blocked reason\" reading of Scaffold Tree viewer should be \"\"", () => readingReads(page, "generate blocked reason", el("Scaffold Tree viewer"), ""));
      await session.step(35, "And the \"generate hint\" reading of Scaffold Tree viewer should be \"Generate from a sample of 750 of 1,000 molecules\"", () => readingReads(page, "generate hint", el("Scaffold Tree viewer"), "Generate from a sample of 750 of 1,000 molecules"));
      await session.step(36, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(37, "Then \"Generate\" icon inside Scaffold Tree viewer should be enabled", () => shouldBe(page, el("\"Generate\" icon inside Scaffold Tree viewer"), "enabled"));
      await session.step(38, "And the \"nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "nodes", el("Scaffold Tree viewer"), 0));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("On 100 molecules the magic wand generates from them all", async () => {
      await session.step(42, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(43, "When user picks \"Chem > Analyze > Scaffold Tree\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Scaffold Tree"));
      await session.step(44, "Then Scaffold Tree viewer should be visible", () => shouldBe(page, el("Scaffold Tree viewer"), "visible"));
      await session.step(45, "And the \"generate blocked reason\" reading of Scaffold Tree viewer should be \"\"", () => readingReads(page, "generate blocked reason", el("Scaffold Tree viewer"), ""));
      await session.step(46, "And the \"generate hint\" reading of Scaffold Tree viewer should be \"Generate from molecular column\"", () => readingReads(page, "generate hint", el("Scaffold Tree viewer"), "Generate from molecular column"));
      await session.step(47, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(48, "Then \"Generate\" icon inside Scaffold Tree viewer should be enabled", () => shouldBe(page, el("\"Generate\" icon inside Scaffold Tree viewer"), "enabled"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With two tables open the menu binds the viewer to the active one (github-3004)", async () => {
      await session.step(52, "Given user opens spgi dataset keeping the first 20 rows as \"tableA\"", () => openDatasetRowsAs(page, ds("spgi"), 20, "tableA"));
      await session.step(53, "And user opens spgi dataset keeping the first 30 rows as \"tableB\"", () => openDatasetRowsAs(page, ds("spgi"), 30, "tableB"));
      await session.step(54, "Then the \"tableB\" view should be current", () => viewIsCurrent(page, "tableB"));
      await session.step(55, "When user picks \"Chem > Analyze > Scaffold Tree\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Scaffold Tree"));
      await session.step(56, "Then Scaffold Tree viewer should be bound to table \"tableB\"", () => boundTable(page, el("Scaffold Tree viewer"), "tableB"));
      await session.step(57, "And the open tableview should have 1 Scaffold Tree viewer", () => viewerCount(page, 1, "Scaffold Tree"));
      await session.step(58, "When user switches to the \"tableA\" table view", () => switchTableView(page, "tableA"));
      await session.step(59, "Then the open tableview should have 0 Scaffold Tree viewers", () => viewerCount(page, 0, "Scaffold Tree"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Coloring a scaffold adds the colors column", async () => {
      await session.step(63, "Given user switches to the \"spgi-100\" table view", () => switchTableView(page, "spgi-100"));
      await session.step(64, "And the stand runs the \"Jupyter\" service", () => standRunsService(page, "Jupyter"));
      await session.step(65, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(66, "And user clicks on \"Generate\" icon inside Scaffold Tree viewer", () => clickOn(page, el("\"Generate\" icon inside Scaffold Tree viewer")));
      await session.step(67, "Then Scaffold Tree viewer should have finished building its tree", () => treeBuilt(page, el("Scaffold Tree viewer")));
      await session.step(68, "And the table should not have a column \"Structure colors\"", () => hasNoColumn(page, "Structure colors"));
      await session.step(69, "When user hovers over the \"node 1\" area of Scaffold Tree viewer", () => hoverArea(page, "node 1", el("Scaffold Tree viewer")));
      await session.step(70, "And user clicks on the \"color icon of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "color icon of node 1", el("Scaffold Tree viewer")));
      await session.step(71, "Then the \"colored nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "colored nodes", el("Scaffold Tree viewer"), 1));
      await session.step(72, "And the table should have a column \"Structure colors\"", () => hasColumn(page, "Structure colors"));
      await session.step(73, "And the \"colors column\" reading of Scaffold Tree viewer should be \"Structure colors\"", () => readingReads(page, "colors column", el("Scaffold Tree viewer"), "Structure colors"));
      await session.step(74, "And \"Structure colors\" column should have at least 1 distinct values", () => distinctValues(page, "Structure colors", 1));
      await session.step(75, "And \"Structure colors\" column should be color-coded categorically", () => colorCodedCategorically(page, "Structure colors"));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A scatter plot colored by the colors column takes the scaffold's colors", async () => {
      await session.step(79, "Given user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(80, "When user picks \"Structure colors\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "Structure colors", "color", el("scatter plot viewer")));
      await session.step(81, "Then \"colorColumnName\" property of scatter plot viewer should be \"Structure colors\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "Structure colors"));
      await session.step(82, "And the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing an uncolored scaffold keeps the colors column and the plot colored by it", async () => {
      await session.step(86, "Given user switches to the \"spgi-100\" table view", () => switchTableView(page, "spgi-100"));
      await session.step(87, "Then the \"colored nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "colored nodes", el("Scaffold Tree viewer"), 1));
      await session.step(88, "When user remembers the \"nodes\" reading of Scaffold Tree viewer", () => rememberReading(page, "nodes", el("Scaffold Tree viewer")));
      await session.step(89, "And user hovers over the \"node 2\" area of Scaffold Tree viewer", () => hoverArea(page, "node 2", el("Scaffold Tree viewer")));
      await session.step(90, "And user clicks on the \"remove icon of node 2\" area of Scaffold Tree viewer", () => clickArea(page, "remove icon of node 2", el("Scaffold Tree viewer")));
      await session.step(91, "And user clicks on \"Yes\" button in \"Remove scaffold\" dialog", () => clickOn(page, el("\"Yes\" button in \"Remove scaffold\" dialog")));
      await session.step(92, "Then the \"nodes\" reading of Scaffold Tree viewer should not be as remembered", () => readingNotAsRemembered(page, "nodes", el("Scaffold Tree viewer")));
      await session.step(93, "And the \"colored nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "colored nodes", el("Scaffold Tree viewer"), 1));
      await session.step(94, "And the table should have a column \"Structure colors\"", () => hasColumn(page, "Structure colors"));
      await session.step(95, "And \"colorColumnName\" property of scatter plot viewer should be \"Structure colors\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "Structure colors"));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the last colored scaffold drops the colors column and the plot's coloring", async () => {
      await session.step(99, "Then the \"color of node 1\" reading of Scaffold Tree viewer should not be \"\"", () => readingDoesNotRead(page, "color of node 1", el("Scaffold Tree viewer"), ""));
      await session.step(100, "When user hovers over the \"node 1\" area of Scaffold Tree viewer", () => hoverArea(page, "node 1", el("Scaffold Tree viewer")));
      await session.step(101, "And user clicks on the \"remove icon of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "remove icon of node 1", el("Scaffold Tree viewer")));
      await session.step(102, "And user clicks on \"Yes\" button in \"Remove scaffold\" dialog", () => clickOn(page, el("\"Yes\" button in \"Remove scaffold\" dialog")));
      await session.step(103, "Then the \"colored nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "colored nodes", el("Scaffold Tree viewer"), 0));
      await session.step(104, "And the table should not have a column \"Structure colors\"", () => hasNoColumn(page, "Structure colors"));
      await session.step(105, "And \"colorColumnName\" property of scatter plot viewer should be \"\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), ""));
      await session.step(106, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
