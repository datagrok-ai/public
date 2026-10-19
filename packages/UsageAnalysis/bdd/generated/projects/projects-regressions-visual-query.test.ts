/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-visual-query.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, GROK-20535]
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
import {addToBuilderRow, builderPreviewRows, builderRowHolds} from '../../bindings/queries.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, isExpanded, pressKey, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clickRunsQuery} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, currentViewType, dialogCloses, doubleClickOpensProject, noProjectOnServer, noQueryOnServer, queriesOnServer, reloadedByDataSync, savedWithDataSync, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a project of a visual query with a parameter", () => {
  const session = feature(test, "features/projects/projects-regressions-visual-query.feature", import.meta.url);
  test("Projects regressions: a project of a visual query with a parameter", {tag: ["@journey", "@realizes:views.projects", "@known-failure:77", "@realizes:GROK-20535"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And no query named \"BDDRegVQ{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDRegVQ{time}")));
    await session.step(24, "And no project named \"BDDRegVQProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegVQProj{time}")));
    await run.scenario("A visual query on entity_types filters name by a parameter", async () => {
      await session.step(27, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(28, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(29, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(30, "And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas tree node inside browse tree")));
      await session.step(31, "And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree")));
      await session.step(32, "When user picks \"New Visual Query...\" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree", () => pickFromContextMenu(page, "New Visual Query...", el("Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree")));
      await session.step(33, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(34, "When user adds \"name\" to the \"Where\" row of the visual query", () => addToBuilderRow(page, "name", "Where"));
      await session.step(35, "And user enters \"Project\" into \"name\" where condition", () => enterInto(page, "Project", el("\"name\" where condition")));
      await session.step(36, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(37, "And user checks \"name\" where parameter checkbox", () => check(page, el("\"name\" where parameter checkbox")));
      await session.step(38, "Then \"name\" where parameter checkbox should be checked", () => shouldBe(page, el("\"name\" where parameter checkbox"), "checked"));
      await session.step(39, "When user adds \"name\" to the \"Order by\" row of the visual query", () => addToBuilderRow(page, "name", "Order by"));
      await session.step(40, "Then the \"Order by\" row of the visual query should hold \"name\"", () => builderRowHolds(page, "Order by", "name"));
      await session.step(41, "And the visual query should preview 1 row", () => builderPreviewRows(page, 1));
      await session.step(42, "When user enters \"BDDRegVQ{time}\" into Name input", () => enterInto(page, session.text("BDDRegVQ{time}"), el("Name input")));
      await session.step(43, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(44, "Then 1 query named \"BDDRegVQ{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDDRegVQ{time}")));
    });
    await run.scenario("Run with another value, the result is saved as a project with Data sync", async () => {
      await session.step(47, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(48, "When user clicks on \"Run query...\" action in toolbox", () => clickOn(page, el("\"Run query...\" action in toolbox")));
      await session.step(49, "Then \"BDDRegVQ{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDRegVQ{time}\" dialog")), "visible"));
      await session.step(50, "And Name text input in \"BDDRegVQ{time}\" dialog should have value \"Project\"", () => shouldHaveValue(page, el(session.text("Name text input in \"BDDRegVQ{time}\" dialog")), "Project"));
      await session.step(51, "When user enters \"Script\" into Name text input in \"BDDRegVQ{time}\" dialog", () => enterInto(page, "Script", el(session.text("Name text input in \"BDDRegVQ{time}\" dialog"))));
      await session.step(52, "And user clicks on OK button in \"BDDRegVQ{time}\" dialog and the query it runs completes", () => clickRunsQuery(page, el(session.text("OK button in \"BDDRegVQ{time}\" dialog"))));
      await session.step(53, "Then the \"BDDRegVQ{time}\" dialog should close", () => dialogCloses(page, session.text("BDDRegVQ{time}")));
      await session.step(54, "And the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(55, "And the table should have 1 row", () => rowCount(page, 1));
      await session.step(56, "And the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
      await session.step(57, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(58, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(59, "And Data sync switch in \"BDDRegVQ{time}\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el(session.text("Data sync switch in \"BDDRegVQ{time}\" project table in \"Save project\" dialog")), "checked"));
      await session.step(60, "When user enters \"BDDRegVQProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegVQProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(61, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(62, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(63, "And \"Share BDDRegVQProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegVQProj{time}\" dialog")), "visible"));
      await session.step(64, "When user clicks on CANCEL button in \"Share BDDRegVQProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegVQProj{time}\" dialog"))));
      await session.step(65, "Then the \"Share BDDRegVQProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegVQProj{time}")));
      await session.step(66, "And the \"BDDRegVQ{time}\" table of the \"BDDRegVQProj{time}\" project should be saved with data sync", () => savedWithDataSync(page, session.text("BDDRegVQ{time}"), session.text("BDDRegVQProj{time}")));
    });
    await run.scenario("The reopened project holds the row of the value it was run with", async () => {
      await session.step(70, "When user closes all views", () => closeAllViews(page));
      await session.step(71, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(72, "And user enters \"BDDRegVQProj{time}\" into gallery search", () => enterInto(page, session.text("BDDRegVQProj{time}"), el("gallery search")));
      await session.step(73, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(74, "And user double-clicks on BDDRegVQProj{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDRegVQProj{time} gallery card"))));
      await session.step(75, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(76, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(77, "And the \"text of cell 1 of name\" reading of grid should be \"Script\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Script"));
    }, {knownFailure: true, line: 77});
    run.finish();
  });
});
