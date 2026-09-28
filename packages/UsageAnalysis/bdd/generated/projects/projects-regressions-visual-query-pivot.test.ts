/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-visual-query-pivot.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, GROK-20273, GROK-20535]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/minimized-viewers.js';
import '../../bindings/projects-copies.js';
import '../../bindings/projects-derived.js';
import '../../bindings/projects-sources.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clearSavedParameters} from '../../bindings/pivot-table.js';
import {noErrorsButPreviewNoise, onlyRowReads} from '../../bindings/projects-regressions.js';
import {addToBuilderRow, builderRowCount, builderRowHolds, cleanQueryLayout, savedQueryFilters} from '../../bindings/queries.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, doubleClickOn, enterInto, isExpanded, pressKey, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, closeAllViews, currentViewType, dialogCloses, noProjectOnServer, noQueryOnServer, queriesOnServer, reloadedByDataSync, savedWithDataSync, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a project of a visual query with a pivot and a parameter", () => {
  const session = feature(test, "features/projects/projects-regressions-visual-query-pivot.feature", import.meta.url);
  test("Projects regressions: a project of a visual query with a pivot and a parameter", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:GROK-20273", "@known-failure", "@realizes:GROK-20535"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(35, "Given user is logged in", () => loggedIn(page));
    await session.step(36, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(37, "And user clears the saved pivot table parameters", () => clearSavedParameters(page));
    await run.scenario("A project of a visual query with a pivot and a parameter reopens without an error", async () => {
      await session.step(41, "Given no query named \"BDDRegVQPivot{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDRegVQPivot{time}")));
      await session.step(42, "And no project named \"BDDRegVQPivotProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegVQPivotProj{time}")));
      await session.step(43, "And the layout saved for the query \"BDDRegVQPivot{time}\" is deleted at the end", () => cleanQueryLayout(page, session.text("BDDRegVQPivot{time}")));
      await session.step(44, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(45, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(46, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(47, "And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas tree node inside browse tree")));
      await session.step(48, "And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree")));
      await session.step(49, "When user picks \"New Visual Query...\" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity-types tree node inside browse tree", () => pickFromContextMenu(page, "New Visual Query...", el("Databases---Postgres---Datagrok---Schemas---public---entity-types tree node inside browse tree")));
      await session.step(50, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(51, "When user adds \"name\" to the \"Where\" row of the visual query", () => addToBuilderRow(page, "name", "Where"));
      await session.step(52, "And user enters \"Project\" into visual query name condition", () => enterInto(page, "Project", el("visual query name condition")));
      await session.step(53, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(54, "And user checks visual query name parameter checkbox", () => check(page, el("visual query name parameter checkbox")));
      await session.step(55, "Then visual query name parameter checkbox should be checked", () => shouldBe(page, el("visual query name parameter checkbox"), "checked"));
      await session.step(56, "When user adds \"name\" to the \"Group by\" row of the visual query", () => addToBuilderRow(page, "name", "Group by"));
      await session.step(57, "And user adds \"id\" to the \"Aggregate\" row of the visual query", () => addToBuilderRow(page, "id", "Aggregate"));
      await session.step(58, "And user adds \"is_package_entity\" to the \"Pivot\" row of the visual query", () => addToBuilderRow(page, "is_package_entity", "Pivot"));
      await session.step(59, "Then the \"Group by\" row of the visual query should hold \"name\"", () => builderRowHolds(page, "Group by", "name"));
      await session.step(60, "And the \"Pivot\" row of the visual query should hold \"is_package_entity\"", () => builderRowHolds(page, "Pivot", "is_package_entity"));
      await session.step(61, "And the visual query should have run to 1 row", () => builderRowCount(page, 1));
      await session.step(62, "When user enters \"BDDRegVQPivot{time}\" into Name input", () => enterInto(page, session.text("BDDRegVQPivot{time}"), el("Name input")));
      await session.step(63, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(64, "Then 1 query named \"BDDRegVQPivot{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDDRegVQPivot{time}")));
      await session.step(65, "And the query \"BDDRegVQPivot{time}\" on the server should filter \"name\" by \"Project\" as a parameter", () => savedQueryFilters(page, session.text("BDDRegVQPivot{time}"), "name", "Project"));
      await session.step(66, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(67, "When user clicks on \"Run query...\" action in toolbox", () => clickOn(page, el("\"Run query...\" action in toolbox")));
      await session.step(68, "Then \"BDDRegVQPivot{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDRegVQPivot{time}\" dialog")), "visible"));
      await session.step(69, "And Name text input in \"BDDRegVQPivot{time}\" dialog should have value \"Project\"", () => shouldHaveValue(page, el(session.text("Name text input in \"BDDRegVQPivot{time}\" dialog")), "Project"));
      await session.step(70, "When user enters \"Script\" into Name text input in \"BDDRegVQPivot{time}\" dialog", () => enterInto(page, "Script", el(session.text("Name text input in \"BDDRegVQPivot{time}\" dialog"))));
      await session.step(71, "And user clicks on OK button in \"BDDRegVQPivot{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BDDRegVQPivot{time}\" dialog"))));
      await session.step(72, "Then the \"BDDRegVQPivot{time}\" dialog should close", () => dialogCloses(page, session.text("BDDRegVQPivot{time}")));
      await session.step(73, "And the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(74, "And the only row of the table should read \"Script\" in the \"name\" column", () => onlyRowReads(page, "Script", "name"));
      await session.step(75, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(76, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(77, "When user enters \"BDDRegVQPivotProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegVQPivotProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(78, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(79, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(80, "And \"Share BDDRegVQPivotProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegVQPivotProj{time}\" dialog")), "visible"));
      await session.step(81, "When user clicks on CANCEL button in \"Share BDDRegVQPivotProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegVQPivotProj{time}\" dialog"))));
      await session.step(82, "Then the \"Share BDDRegVQPivotProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegVQPivotProj{time}")));
      await session.step(83, "And the \"BDDRegVQPivot{time}\" table of the \"BDDRegVQPivotProj{time}\" project should be saved with data sync", () => savedWithDataSync(page, session.text("BDDRegVQPivot{time}"), session.text("BDDRegVQPivotProj{time}")));
      await session.step(84, "And no errors but the project preview's should have been logged", () => noErrorsButPreviewNoise(page));
      await session.step(85, "When user closes all views", () => closeAllViews(page));
      await session.step(86, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(87, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(88, "And user enters \"BDDRegVQPivotProj{time}\" into gallery search", () => enterInto(page, session.text("BDDRegVQPivotProj{time}"), el("gallery search")));
      await session.step(89, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(90, "Given user watches the task bar", () => watchTaskBar(page));
      await session.step(91, "When user double-clicks on BDDRegVQPivotProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDRegVQPivotProj{time} gallery card"))));
      await session.step(92, "Then the task bar should have finished \"Opening project\"", () => taskBarFinished(page, "Opening project"));
      await session.step(93, "And the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(94, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(95, "And the table should have 1 row", () => rowCount(page, 1));
      await session.step(96, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The reopened pivoted visual query holds the row of the value it was run with", async () => {
      await session.step(101, "Then the only row of the table should read \"Script\" in the \"name\" column", () => onlyRowReads(page, "Script", "name"));
    }, {knownFailure: true});
    run.finish();
  });
});
