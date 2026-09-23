/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scaffold-tree/colors-and-limits.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.scaffold-tree-add-filter]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {treeBuilt} from '../../bindings/scaffold-tree.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {distinctValues, hasColumn, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCodedCategorically} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openDatasetRowsAs, switchTableView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, clickArea, hoverArea, legendLists, noErrors, painted, propertyShouldBe, readingAtLeast, readingIs, readingLower, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector, readingIncludes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scaffold Tree colors, blocked generation and two tables", () => {
  const session = feature(test, "features/scaffold-tree/colors-and-limits.feature", import.meta.url);
  test("Scaffold Tree colors, blocked generation and two tables", {tag: ["@journey", "@realizes:chem.cp.scaffold-tree-add-filter", "@known-failure", "@GROK-18286"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(14, "And user picks \"Chem > Analyze > Scaffold Tree\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Scaffold Tree"));
    await run.scenario("Coloring a scaffold adds the colors column", async () => {
      await session.step(17, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(18, "And user clicks on \"Generate\" icon inside Scaffold Tree viewer", () => clickOn(page, el("\"Generate\" icon inside Scaffold Tree viewer")));
      await session.step(19, "Then Scaffold Tree viewer should have finished building its tree", () => treeBuilt(page, el("Scaffold Tree viewer")));
      await session.step(20, "And the table should not have a column \"Structure colors\"", () => hasNoColumn(page, "Structure colors"));
      await session.step(21, "When user hovers over the \"node 1\" area of Scaffold Tree viewer", () => hoverArea(page, "node 1", el("Scaffold Tree viewer")));
      await session.step(22, "And user clicks on the \"color icon of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "color icon of node 1", el("Scaffold Tree viewer")));
      await session.step(23, "Then the \"colored nodes\" reading of Scaffold Tree viewer should be 1", () => readingIs(page, "colored nodes", el("Scaffold Tree viewer"), 1));
      await session.step(24, "And the table should have a column \"Structure colors\"", () => hasColumn(page, "Structure colors"));
      await session.step(25, "And the \"colors column\" reading of Scaffold Tree viewer should be \"Structure colors\"", () => readingReads(page, "colors column", el("Scaffold Tree viewer"), "Structure colors"));
      await session.step(26, "And \"Structure colors\" column should have at least 1 distinct values", () => distinctValues(page, "Structure colors", 1));
      await session.step(27, "And \"Structure colors\" column should be color-coded categorically", () => colorCodedCategorically(page, "Structure colors"));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A scatter plot colored by the colors column takes the scaffold's colors", async () => {
      await session.step(31, "Given user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(32, "When user picks \"Structure colors\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "Structure colors", "color", el("scatter plot viewer")));
      await session.step(33, "Then \"colorColumnName\" property of scatter plot viewer should be \"Structure colors\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "Structure colors"));
      await session.step(34, "And the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The magic wand is blocked on a table with no molecule column", async () => {
      await session.step(38, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
      await session.step(39, "And user adds Scaffold Tree viewer", () => addViewer(page, "Scaffold Tree"));
      await session.step(40, "Then Scaffold Tree viewer should be visible", () => shouldBe(page, el("Scaffold Tree viewer"), "visible"));
      await session.step(41, "And the \"generate blocked reason\" reading of Scaffold Tree viewer should be \"There is no molecule column in the table\"", () => readingReads(page, "generate blocked reason", el("Scaffold Tree viewer"), "There is no molecule column in the table"));
      await session.step(42, "And the \"message\" reading of Scaffold Tree viewer should include the text \"No molecule column found\"", () => readingIncludes(page, "message", el("Scaffold Tree viewer"), "No molecule column found"));
      await session.step(43, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(44, "Then \"Generate\" icon inside Scaffold Tree viewer should be disabled", () => shouldBe(page, el("\"Generate\" icon inside Scaffold Tree viewer"), "disabled"));
      await session.step(45, "And the \"nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "nodes", el("Scaffold Tree viewer"), 0));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The magic wand is blocked past 500 structure categories", async () => {
      await session.step(49, "Given user opens mol1K dataset", () => openDataset(page, ds("mol1K")));
      await session.step(50, "When user picks \"Chem > Analyze > Scaffold Tree\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Scaffold Tree"));
      await session.step(51, "Then Scaffold Tree viewer should be visible", () => shouldBe(page, el("Scaffold Tree viewer"), "visible"));
      await session.step(52, "And the \"generate blocked reason\" reading of Scaffold Tree viewer should be \"The number of molecules exceeds the limit of 500\"", () => readingReads(page, "generate blocked reason", el("Scaffold Tree viewer"), "The number of molecules exceeds the limit of 500"));
      await session.step(53, "When user hovers over Scaffold Tree viewer", () => hoverOver(page, el("Scaffold Tree viewer")));
      await session.step(54, "Then \"Generate\" icon inside Scaffold Tree viewer should be disabled", () => shouldBe(page, el("\"Generate\" icon inside Scaffold Tree viewer"), "disabled"));
      await session.step(55, "And the \"nodes\" reading of Scaffold Tree viewer should be 0", () => readingIs(page, "nodes", el("Scaffold Tree viewer"), 0));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With two tables open the menu binds the viewer to the active one (github-3004)", async () => {
      await session.step(59, "Given user opens spgi dataset keeping the first 20 rows as \"tableA\"", () => openDatasetRowsAs(page, ds("spgi"), 20, "tableA"));
      await session.step(60, "And user opens spgi dataset keeping the first 30 rows as \"tableB\"", () => openDatasetRowsAs(page, ds("spgi"), 30, "tableB"));
      await session.step(61, "Then the \"tableB\" view should be current", () => viewIsCurrent(page, "tableB"));
      await session.step(62, "When user picks \"Chem > Analyze > Scaffold Tree\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Scaffold Tree"));
      await session.step(63, "Then Scaffold Tree viewer should be bound to table \"tableB\"", () => boundTable(page, el("Scaffold Tree viewer"), "tableB"));
      await session.step(64, "And the open tableview should have 1 Scaffold Tree viewer", () => viewerCount(page, 1, "Scaffold Tree"));
      await session.step(65, "When user switches to the \"tableA\" table view", () => switchTableView(page, "tableA"));
      await session.step(66, "Then the open tableview should have 0 Scaffold Tree viewers", () => viewerCount(page, 0, "Scaffold Tree"));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing a colored scaffold leaves the scatter plot that is colored by it", async () => {
      await session.step(71, "Given user switches to the \"spgi-100\" table view", () => switchTableView(page, "spgi-100"));
      await session.step(72, "Then the \"nodes\" reading of Scaffold Tree viewer should be at least 1", () => readingAtLeast(page, "nodes", el("Scaffold Tree viewer"), 1));
      await session.step(73, "When user hovers over the \"node 1\" area of Scaffold Tree viewer", () => hoverArea(page, "node 1", el("Scaffold Tree viewer")));
      await session.step(74, "And user clicks on the \"remove icon of node 1\" area of Scaffold Tree viewer", () => clickArea(page, "remove icon of node 1", el("Scaffold Tree viewer")));
      await session.step(75, "And user clicks on \"Yes\" button in \"Remove scaffold\" dialog", () => clickOn(page, el("\"Yes\" button in \"Remove scaffold\" dialog")));
      await session.step(76, "Then the \"nodes\" reading of Scaffold Tree viewer should be lower than before", () => readingLower(page, "nodes", el("Scaffold Tree viewer")));
      await session.step(77, "And \"colorColumnName\" property of scatter plot viewer should be \"Structure colors\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "Structure colors"));
      await session.step(78, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true});
    run.finish();
  });
});
