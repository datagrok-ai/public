/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/parameterized-queries.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.queries]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, followingShouldBe, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, currentViewType, dialogCloses, standHasReachableConnection, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Running a query with parameters", () => {
  const session = feature(test, "features/queries/parameterized-queries.feature", import.meta.url);
  test("The Orders query asks for its eight typed parameters", {tag: ["@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the stand has a reachable \"PostgresTest\" connection", () => standHasReachableConnection(page, "PostgresTest"));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(24, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(25, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(26, "When user picks \"Run\" from the context menu of Databases---Postgres---NorthwindTest---Orders tree node inside browse tree", () => pickFromContextMenu(page, "Run", el("Databases---Postgres---NorthwindTest---Orders tree node inside browse tree")));
    await session.step(27, "Then \"Orders\" dialog should be visible", () => shouldBe(page, el("\"Orders\" dialog"), "visible"));
    await session.step(28, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Employee Id\" input in \"Orders\" dialog"],["\"Ship Via\" input in \"Orders\" dialog"],["Freight input in \"Orders\" dialog"],["\"Ship Country\" input in \"Orders\" dialog"],["\"Ship City\" input in \"Orders\" dialog"],["\"Freight Less1000\" input in \"Orders\" dialog"],["\"Required Date\" input in \"Orders\" dialog"],["\"Order Date\" input in \"Orders\" dialog"]]), [["\"Employee Id\" input in \"Orders\" dialog"],["\"Ship Via\" input in \"Orders\" dialog"],["Freight input in \"Orders\" dialog"],["\"Ship Country\" input in \"Orders\" dialog"],["\"Ship City\" input in \"Orders\" dialog"],["\"Freight Less1000\" input in \"Orders\" dialog"],["\"Required Date\" input in \"Orders\" dialog"],["\"Order Date\" input in \"Orders\" dialog"]]);
    await session.step(37, "When user clicks on CANCEL button in \"Orders\" dialog", () => clickOn(page, el("CANCEL button in \"Orders\" dialog")));
    await session.step(38, "Then the \"Orders\" dialog should close", () => dialogCloses(page, "Orders"));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
    await session.step(40, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A string choice runs the query, and REFRESH runs it again with another", {tag: ["@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the stand has a reachable \"PostgresTest\" connection", () => standHasReachableConnection(page, "PostgresTest"));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(43, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(44, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(45, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(46, "When user picks \"Run\" from the context menu of Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree", () => pickFromContextMenu(page, "Run", el("Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree")));
    await session.step(47, "Then \"PostgresByStringChoices\" dialog should be visible", () => shouldBe(page, el("\"PostgresByStringChoices\" dialog"), "visible"));
    await session.step(48, "When user selects \"France\" in \"Ship Country\" input in \"PostgresByStringChoices\" dialog", () => selectIn(page, "France", el("\"Ship Country\" input in \"PostgresByStringChoices\" dialog")));
    await session.step(49, "And user clicks on OK button in \"PostgresByStringChoices\" dialog", () => clickOn(page, el("OK button in \"PostgresByStringChoices\" dialog")));
    await session.step(50, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
    await session.step(51, "And the table should have 77 rows", () => rowCount(page, 77));
    await session.step(52, "And every value of \"shipcountry\" column should match \"^France$\"", () => everyValueMatches(page, "shipcountry", "^France$"));
    await session.step(53, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(54, "When user selects \"USA\" in \"Ship Country\" input in toolbox", () => selectIn(page, "USA", el("\"Ship Country\" input in toolbox")));
    await session.step(55, "And user clicks on REFRESH button in toolbox", () => clickOn(page, el("REFRESH button in toolbox")));
    await session.step(56, "Then the table should have 122 rows", () => rowCount(page, 122));
    await session.step(57, "And no errors should have been logged", () => noErrors(page));
    await session.step(58, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
