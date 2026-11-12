/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/parameterized-queries.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.queries, core#641]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {appendToEditor, clickOn, enterInto, followingShouldBe, isExpanded, replaceCode, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, currentViewType, dialogCloses, noQueryOnServer, standHasReachableConnection, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Running a query with parameters", () => {
  const session = feature(test, "features/queries/parameterized-queries.feature", import.meta.url);
  test("The Orders query asks for its eight typed parameters", {tag: ["@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "Given the stand has a reachable \"PostgresTest\" connection", () => standHasReachableConnection(page, "PostgresTest"));
    await session.step(32, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(33, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(34, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(35, "When user picks \"Run\" from the context menu of Databases---Postgres---NorthwindTest---Orders tree node inside browse tree", () => pickFromContextMenu(page, "Run", el("Databases---Postgres---NorthwindTest---Orders tree node inside browse tree")));
    await session.step(36, "Then \"Orders\" dialog should be visible", () => shouldBe(page, el("\"Orders\" dialog"), "visible"));
    await session.step(37, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Employee Id\" input in \"Orders\" dialog"],["\"Ship Via\" input in \"Orders\" dialog"],["Freight input in \"Orders\" dialog"],["\"Ship Country\" input in \"Orders\" dialog"],["\"Ship City\" input in \"Orders\" dialog"],["\"Freight Less1000\" input in \"Orders\" dialog"],["\"Required Date\" input in \"Orders\" dialog"],["\"Order Date\" input in \"Orders\" dialog"]]), [["\"Employee Id\" input in \"Orders\" dialog"],["\"Ship Via\" input in \"Orders\" dialog"],["Freight input in \"Orders\" dialog"],["\"Ship Country\" input in \"Orders\" dialog"],["\"Ship City\" input in \"Orders\" dialog"],["\"Freight Less1000\" input in \"Orders\" dialog"],["\"Required Date\" input in \"Orders\" dialog"],["\"Order Date\" input in \"Orders\" dialog"]]);
    await session.step(46, "When user clicks on CANCEL button in \"Orders\" dialog", () => clickOn(page, el("CANCEL button in \"Orders\" dialog")));
    await session.step(47, "Then the \"Orders\" dialog should close", () => dialogCloses(page, "Orders"));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
    await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A string choice runs the query, and REFRESH runs it again with another", {tag: ["@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(52, "Given the stand has a reachable \"PostgresTest\" connection", () => standHasReachableConnection(page, "PostgresTest"));
    await session.step(53, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(54, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(55, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(56, "When user picks \"Run\" from the context menu of Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree", () => pickFromContextMenu(page, "Run", el("Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree")));
    await session.step(57, "Then \"PostgresByStringChoices\" dialog should be visible", () => shouldBe(page, el("\"PostgresByStringChoices\" dialog"), "visible"));
    await session.step(58, "When user selects \"France\" in \"Ship Country\" input in \"PostgresByStringChoices\" dialog", () => selectIn(page, "France", el("\"Ship Country\" input in \"PostgresByStringChoices\" dialog")));
    await session.step(59, "And user clicks on OK button in \"PostgresByStringChoices\" dialog", () => clickOn(page, el("OK button in \"PostgresByStringChoices\" dialog")));
    await session.step(60, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
    await session.step(61, "And the table should have 77 rows", () => rowCount(page, 77));
    await session.step(62, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(63, "When user selects \"USA\" in \"Ship Country\" input in toolbox", () => selectIn(page, "USA", el("\"Ship Country\" input in toolbox")));
    await session.step(64, "And user clicks on REFRESH button in toolbox", () => clickOn(page, el("REFRESH button in toolbox")));
    await session.step(65, "Then the table should have 122 rows", () => rowCount(page, 122));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
    await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Run query... on an unsaved parameterized query asks for its parameters", {tag: ["@realizes:views.queries", "@realizes:core#641"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(71, "Given no query named \"BDDQUnsaved{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDQUnsaved{time}")));
    await session.step(72, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(73, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(74, "When user picks \"New Query...\" from the context menu of Databases---Postgres---Datagrok tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(75, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
    await session.step(76, "When user enters \"BDDQUnsaved{time}\" into Name input", () => enterInto(page, session.text("BDDQUnsaved{time}"), el("Name input")));
    await session.step(77, "And user replaces the code of code editor with '--input: string typeName = \"Project\"'", () => replaceCode(page, el("code editor"), "--input: string typeName = \"Project\""));
    await session.step(78, "And user appends \"select name from entity_types where name = @typeName\" to code editor", () => appendToEditor(page, "select name from entity_types where name = @typeName", el("code editor")));
    await session.step(79, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(80, "When user clicks on \"Run query...\" action in toolbox", () => clickOn(page, el("\"Run query...\" action in toolbox")));
    await session.step(81, "Then \"BDDQUnsaved{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDQUnsaved{time}\" dialog")), "visible"));
    await session.step(82, "When user clicks on OK button in \"BDDQUnsaved{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BDDQUnsaved{time}\" dialog"))));
    await session.step(83, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
    await session.step(84, "And the table should have 1 row", () => rowCount(page, 1));
    await session.step(85, "And no errors should have been logged", () => noErrors(page));
  });
});
