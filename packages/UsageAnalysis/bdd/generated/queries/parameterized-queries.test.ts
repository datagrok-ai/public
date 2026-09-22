/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/parameterized-queries.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.queries]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, followingShouldBe, hoverOver, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, currentViewType, dialogCloses, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Running a query with parameters", () => {
  const session = feature(test, "features/queries/parameterized-queries.feature", import.meta.url);
  test("The Orders query asks for its eight typed parameters", {tag: ["@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(25, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(26, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(28, "When user hovers over Databases---Postgres---NorthwindTest---Orders tree node inside browse tree", () => hoverOver(page, el("Databases---Postgres---NorthwindTest---Orders tree node inside browse tree")));
    await session.step(29, "And user picks \"Run\" from the context menu of Databases---Postgres---NorthwindTest---Orders tree node inside browse tree", () => pickFromContextMenu(page, "Run", el("Databases---Postgres---NorthwindTest---Orders tree node inside browse tree")));
    await session.step(30, "Then \"Orders\" dialog should be visible", () => shouldBe(page, el("\"Orders\" dialog"), "visible"));
    await session.step(31, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Employee Id\" input in \"Orders\" dialog"],["\"Ship Via\" input in \"Orders\" dialog"],["Freight input in \"Orders\" dialog"],["\"Ship Country\" input in \"Orders\" dialog"],["\"Ship City\" input in \"Orders\" dialog"],["\"Freight Less1000\" input in \"Orders\" dialog"],["\"Required Date\" input in \"Orders\" dialog"],["\"Order Date\" input in \"Orders\" dialog"]]), [["\"Employee Id\" input in \"Orders\" dialog"],["\"Ship Via\" input in \"Orders\" dialog"],["Freight input in \"Orders\" dialog"],["\"Ship Country\" input in \"Orders\" dialog"],["\"Ship City\" input in \"Orders\" dialog"],["\"Freight Less1000\" input in \"Orders\" dialog"],["\"Required Date\" input in \"Orders\" dialog"],["\"Order Date\" input in \"Orders\" dialog"]]);
    await session.step(40, "When user clicks on CANCEL button in \"Orders\" dialog", () => clickOn(page, el("CANCEL button in \"Orders\" dialog")));
    await session.step(41, "Then the \"Orders\" dialog should close", () => dialogCloses(page, "Orders"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
    await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A string choice runs the query, and REFRESH runs it again with another", {tag: ["@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(46, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(47, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(48, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(50, "When user hovers over Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree", () => hoverOver(page, el("Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree")));
    await session.step(51, "And user picks \"Run\" from the context menu of Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree", () => pickFromContextMenu(page, "Run", el("Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree")));
    await session.step(52, "Then \"PostgresByStringChoices\" dialog should be visible", () => shouldBe(page, el("\"PostgresByStringChoices\" dialog"), "visible"));
    await session.step(53, "When user selects \"France\" in \"Ship Country\" input in \"PostgresByStringChoices\" dialog", () => selectIn(page, "France", el("\"Ship Country\" input in \"PostgresByStringChoices\" dialog")));
    await session.step(54, "And user clicks on OK button in \"PostgresByStringChoices\" dialog", () => clickOn(page, el("OK button in \"PostgresByStringChoices\" dialog")));
    await session.step(55, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
    await session.step(56, "And the table should have 77 rows", () => rowCount(page, 77));
    await session.step(57, "And every value of \"shipcountry\" column should match \"^France$\"", () => everyValueMatches(page, "shipcountry", "^France$"));
    await session.step(58, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(59, "When user selects \"USA\" in \"Ship Country\" input in toolbox", () => selectIn(page, "USA", el("\"Ship Country\" input in toolbox")));
    await session.step(60, "And user clicks on REFRESH button in toolbox", () => clickOn(page, el("REFRESH button in toolbox")));
    await session.step(61, "Then the table should have 122 rows", () => rowCount(page, 122));
    await session.step(62, "And no errors should have been logged", () => noErrors(page));
    await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("CHEMBL's FRAC search asks for a mechanism and a substructure", {tag: ["@realizes:views.queries", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(67, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(68, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(69, "And Databases---Postgres---CHEMBL tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---CHEMBL tree node inside browse tree")));
    await session.step(71, "And \"Databases---Postgres---CHEMBL---Search-\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"Databases---Postgres---CHEMBL---Search-\" tree node inside browse tree")));
    await session.step(73, "When user hovers over Databases---Postgres---CHEMBL---Search-----By-FRAC-Classification-And-Substructure tree node inside browse tree", () => hoverOver(page, el("Databases---Postgres---CHEMBL---Search-----By-FRAC-Classification-And-Substructure tree node inside browse tree")));
    await session.step(74, "And user picks \"Run\" from the context menu of Databases---Postgres---CHEMBL---Search-----By-FRAC-Classification-And-Substructure tree node inside browse tree", () => pickFromContextMenu(page, "Run", el("Databases---Postgres---CHEMBL---Search-----By-FRAC-Classification-And-Substructure tree node inside browse tree")));
    await session.step(75, "Then \"Search | By FRAC Classification And Substructure\" dialog should be visible", () => shouldBe(page, el("\"Search | By FRAC Classification And Substructure\" dialog"), "visible"));
    await session.step(76, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Mechanism input in \"Search | By FRAC Classification And Substructure\" dialog"],["Substructure input in \"Search | By FRAC Classification And Substructure\" dialog"]]), [["Mechanism input in \"Search | By FRAC Classification And Substructure\" dialog"],["Substructure input in \"Search | By FRAC Classification And Substructure\" dialog"]]);
    await session.step(79, "When user clicks on OK button in \"Search | By FRAC Classification And Substructure\" dialog", () => clickOn(page, el("OK button in \"Search | By FRAC Classification And Substructure\" dialog")));
    await session.step(80, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
    await session.step(81, "And the table should have 26 rows", () => rowCount(page, 26));
    await session.step(82, "And no errors should have been logged", () => noErrors(page));
    await session.step(83, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
