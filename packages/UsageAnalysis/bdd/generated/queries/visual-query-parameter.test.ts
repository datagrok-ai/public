/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/visual-query-parameter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.queries, sharing.share-dialog]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {addToBuilderRow, builderPreviewRows, builderRowHidden, builderRowHolds, builderRowOffers, builderRowShown, clickBuilderPlus, pickJoinedColumn, removeFromBuilderRow} from '../../bindings/queries.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {appendToEditor, check, clickOn, enterInto, isExpanded, shouldBe, shouldContainText, shouldHaveValue, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {clickRunsQuery} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, currentViewType, dialogCloses, entityHasLayout, layoutsDeleted, noQueryOnServer, pickSharingUser, queriesOnServer, queryPostProcess, refreshBrowse, sharingUserLevel, toolboxPaneHidden, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, infoBalloonText, menuLists, noBalloons, noErrors, pickFromContextMenu, pickFromOpenMenu, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A visual query with a parameter, from the builder to the shared query", () => {
  const session = feature(test, "features/queries/visual-query-parameter.feature", import.meta.url);
  test("A visual query with a parameter, from the builder to the shared query", {tag: ["@journey", "@serial", "@realizes:views.queries", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "And no query named \"BDD-VQ-{time}\" is on the server", () => noQueryOnServer(page, session.text("BDD-VQ-{time}")));
    await session.step(30, "And the layouts named \"BDD-VQ-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("BDD-VQ-{time}")));
    await run.scenario("The Where condition is exposed as a function parameter", async () => {
      await session.step(33, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(34, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(35, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(36, "And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas tree node inside browse tree")));
      await session.step(37, "And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree")));
      await session.step(38, "When user picks \"New Visual Query...\" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree", () => pickFromContextMenu(page, "New Visual Query...", el("Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree")));
      await session.step(39, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(40, "When user adds \"name\" to the \"Where\" row of the visual query", () => addToBuilderRow(page, "name", "Where"));
      await session.step(41, "Then the \"Where\" row of the visual query should hold \"name\"", () => builderRowHolds(page, "Where", "name"));
      await session.step(42, "When user enters \"User\" into \"name\" where condition", () => enterInto(page, "User", el("\"name\" where condition")));
      await session.step(43, "Then the visual query should preview 1 row", () => builderPreviewRows(page, 1));
      await session.step(44, "When user checks \"name\" where parameter checkbox", () => check(page, el("\"name\" where parameter checkbox")));
      await session.step(45, "Then no errors should have been logged", () => noErrors(page));
      await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("An aggregation brings up Pivot and Having, and Order by offers only the groups", async () => {
      await session.step(49, "Then the \"Pivot\" row of the visual query should be hidden", () => builderRowHidden(page, "Pivot"));
      await session.step(50, "And the \"Having\" row of the visual query should be hidden", () => builderRowHidden(page, "Having"));
      await session.step(51, "When user adds \"name\" to the \"Group by\" row of the visual query", () => addToBuilderRow(page, "name", "Group by"));
      await session.step(52, "And user adds \"id\" to the \"Aggregate\" row of the visual query", () => addToBuilderRow(page, "id", "Aggregate"));
      await session.step(53, "Then the \"Aggregate\" row of the visual query should hold \"values(id)\"", () => builderRowHolds(page, "Aggregate", "values(id)"));
      await session.step(54, "And the \"Pivot\" row of the visual query should be shown", () => builderRowShown(page, "Pivot"));
      await session.step(55, "And the \"Having\" row of the visual query should be shown", () => builderRowShown(page, "Having"));
      await session.step(56, "And the \"Order by\" row of the visual query should offer the columns \"name\"", () => builderRowOffers(page, "Order by", "name"));
      await session.step(57, "When user adds \"is_package_entity\" to the \"Pivot\" row of the visual query", () => addToBuilderRow(page, "is_package_entity", "Pivot"));
      await session.step(58, "Then the \"Pivot\" row of the visual query should hold \"is_package_entity\"", () => builderRowHolds(page, "Pivot", "is_package_entity"));
      await session.step(59, "When user adds \"id\" to the \"Having\" row of the visual query", () => addToBuilderRow(page, "id", "Having"));
      await session.step(60, "And user enters \"> 1\" into \"values-id\" having condition", () => enterInto(page, "> 1", el("\"values-id\" having condition")));
      await session.step(61, "Then the visual query should preview 0 rows", () => builderPreviewRows(page, 0));
      await session.step(62, "When user enters \"1\" into \"values-id\" having condition", () => enterInto(page, "1", el("\"values-id\" having condition")));
      await session.step(63, "Then the visual query should preview 1 row", () => builderPreviewRows(page, 1));
      await session.step(64, "When user adds \"name\" to the \"Order by\" row of the visual query", () => addToBuilderRow(page, "name", "Order by"));
      await session.step(65, "Then the \"Order by\" row of the visual query should hold \"name\"", () => builderRowHolds(page, "Order by", "name"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The saved query asks for its parameter, and REFRESH runs it with another", async () => {
      await session.step(70, "When user enters \"BDD-VQ-{time}\" into Name input", () => enterInto(page, session.text("BDD-VQ-{time}"), el("Name input")));
      await session.step(71, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(72, "Then 1 query named \"BDD-VQ-{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDD-VQ-{time}")));
      await session.step(73, "When user closes all views", () => closeAllViews(page));
      await session.step(74, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(75, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(76, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(77, "And user picks \"Run\" from the context menu of Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Run", el(session.text("Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree"))));
      await session.step(78, "Then \"BDD-VQ-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-VQ-{time}\" dialog")), "visible"));
      await session.step(79, "And Name input in \"BDD-VQ-{time}\" dialog should have value \"User\"", () => shouldHaveValue(page, el(session.text("Name input in \"BDD-VQ-{time}\" dialog")), "User"));
      await session.step(80, "When user clicks on OK button in \"BDD-VQ-{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BDD-VQ-{time}\" dialog"))));
      await session.step(81, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(82, "And the table should have 1 row", () => rowCount(page, 1));
      await session.step(83, "And the table should have the columns \"name, false\"", () => columnsExactly(page, "name, false"));
      await session.step(84, "And the value of \"name\" column in row 1 should be \"User\"", () => valueInRow(page, "name", 1, "User"));
      await session.step(85, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(86, "When user types \"UserGroup\" into Name input in toolbox", () => typeInto(page, "UserGroup", el("Name input in toolbox")));
      await session.step(87, "And user clicks on REFRESH button in toolbox and the query it runs completes", () => clickRunsQuery(page, el("REFRESH button in toolbox")));
      await session.step(88, "Then the value of \"name\" column in row 1 should be \"UserGroup\"", () => valueInRow(page, "name", 1, "UserGroup"));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
      await session.step(90, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Edit changes the query, and adds a post-process and a layout", async () => {
      await session.step(93, "When user closes all views", () => closeAllViews(page));
      await session.step(94, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
      await session.step(95, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(96, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(97, "When user picks \"Edit...\" from the context menu of Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree"))));
      await session.step(98, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(99, "And the \"Where\" row of the visual query should hold \"name\"", () => builderRowHolds(page, "Where", "name"));
      await session.step(100, "And \"name\" where condition should have value \"User\"", () => shouldHaveValue(page, el("\"name\" where condition"), "User"));
      await session.step(101, "And \"name\" where parameter checkbox should be checked", () => shouldBe(page, el("\"name\" where parameter checkbox"), "checked"));
      await session.step(102, "And the \"Pivot\" row of the visual query should hold \"is_package_entity\"", () => builderRowHolds(page, "Pivot", "is_package_entity"));
      await session.step(103, "When user removes \"is_package_entity\" from the \"Pivot\" row of the visual query", () => removeFromBuilderRow(page, "is_package_entity", "Pivot"));
      await session.step(104, "Then the \"Pivot\" row of the visual query should hold \"\"", () => builderRowHolds(page, "Pivot", ""));
      await session.step(105, "When user clicks on Post-Process tab", () => clickOn(page, el("Post-Process tab")));
      await session.step(106, "And user appends \"grok.shell.info('VQ' + result.rowCount);\" to code editor", () => appendToEditor(page, "grok.shell.info('VQ' + result.rowCount);", el("code editor")));
      await session.step(109, "And user clicks on Layout tab", () => clickOn(page, el("Layout tab")));
      await session.step(110, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(111, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(112, "When user clicks on bar chart icon in toolbox", () => clickOn(page, el("bar chart icon in toolbox")));
      await session.step(113, "Then bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
      await session.step(114, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(115, "Then the query \"BDD-VQ-{time}\" on the server should have a layout", () => entityHasLayout(page, "query", session.text("BDD-VQ-{time}")));
      await session.step(116, "And the query \"BDD-VQ-{time}\" on the server should have a post-process containing \"grok.shell.info('VQ' + result.rowCount);\"", () => queryPostProcess(page, session.text("BDD-VQ-{time}"), "grok.shell.info('VQ' + result.rowCount);"));
      await session.step(117, "And no errors should have been logged", () => noErrors(page));
      await session.step(118, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The owner shares the visual query", async () => {
      await session.step(121, "When user closes all views", () => closeAllViews(page));
      await session.step(122, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
      await session.step(123, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(124, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(125, "When user picks \"Share...\" from the context menu of Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el(session.text("Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree"))));
      await session.step(126, "Then \"Share BDD-VQ-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD-VQ-{time}\" dialog")), "visible"));
      await session.step(127, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(128, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDD-VQ-{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDD-VQ-{time}\" dialog"))));
      await session.step(129, "And user unchecks \"Send notifications\" input in \"Share BDD-VQ-{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDD-VQ-{time}\" dialog"))));
      await session.step(130, "And user clicks on OK button in \"Share BDD-VQ-{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDD-VQ-{time}\" dialog"))));
      await session.step(131, "Then the \"Share BDD-VQ-{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDD-VQ-{time}")));
      await session.step(132, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(133, "When user picks \"Share...\" from the context menu of Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el(session.text("Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree"))));
      await session.step(134, "Then \"Share BDD-VQ-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD-VQ-{time}\" dialog")), "visible"));
      await session.step(135, "And the access level of the sharing user in \"Share BDD-VQ-{time}\" dialog should be \"View and use\"", () => sharingUserLevel(page, el(session.text("\"Share BDD-VQ-{time}\" dialog")), "View and use"));
      await session.step(136, "When user clicks on CANCEL button in \"Share BDD-VQ-{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDD-VQ-{time}\" dialog"))));
      await session.step(137, "Then the \"Share BDD-VQ-{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDD-VQ-{time}")));
      await session.step(138, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A click on the saved query previews the edited query with its layout and runs its post-process", async () => {
      await session.step(142, "When user clicks on Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree", () => clickOn(page, el(session.text("Databases---Postgres---Datagrok---BDD-VQ-{time} tree node inside browse tree"))));
      await session.step(143, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(144, "And the table should have the columns \"name, values(id)\"", () => columnsExactly(page, "name, values(id)"));
      await session.step(145, "And the open tableview should have 1 bar chart viewer", () => viewerCount(page, 1, "bar chart"));
      await session.step(146, "And an info balloon containing \"VQ1\" should have been shown", () => infoBalloonText(page, "VQ1"));
      await session.step(147, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With a table joined, the + of a row offers the columns table by table", async () => {
      await session.step(150, "When user closes all views", () => closeAllViews(page));
      await session.step(151, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
      await session.step(152, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(153, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(154, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(155, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(156, "And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas tree node inside browse tree")));
      await session.step(157, "And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree")));
      await session.step(158, "When user picks \"New Visual Query...\" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree", () => pickFromContextMenu(page, "New Visual Query...", el("Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree")));
      await session.step(159, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(160, "When user clicks on the + of the \"Data\" row of the visual query", () => clickBuilderPlus(page, "Data"));
      await session.step(161, "And user picks \"datagrok > public > entity_types_permissions\" from the open menu", () => pickFromOpenMenu(page, "datagrok > public > entity_types_permissions"));
      await session.step(162, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(163, "And user clicks on the + of the \"Group by\" row of the visual query", () => clickBuilderPlus(page, "Group by"));
      await session.step(164, "Then the open menu should list \"datagrok.public.entity_types\"", () => menuLists(page, "datagrok.public.entity_types"));
      await session.step(165, "And the open menu should list \"datagrok.public.entity_types_permissions\"", () => menuLists(page, "datagrok.public.entity_types_permissions"));
      await session.step(166, "When user picks the \"entity_type_id\" column of the \"datagrok.public.entity_types_permissions\" table from the open menu", () => pickJoinedColumn(page, "entity_type_id", "datagrok.public.entity_types_permissions"));
      await session.step(167, "Then the \"Group by\" row of the visual query should hold \"datagrok.public.entity_types_permissions.entity_type_id\"", () => builderRowHolds(page, "Group by", "datagrok.public.entity_types_permissions.entity_type_id"));
      await session.step(168, "And no errors should have been logged", () => noErrors(page));
      await session.step(169, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(170, "When user closes all views", () => closeAllViews(page));
    });
    run.finish();
  });
});
