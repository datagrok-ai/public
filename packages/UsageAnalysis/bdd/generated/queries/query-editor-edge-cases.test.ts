/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/query-editor-edge-cases.feature
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
import {appendToEditor, clickOn, enterInto, isExpanded, replaceCode, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, currentViewType, dialogCloses, noQueryOnServer, queriesOnServer, queryOnServer, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Query editor and Rename edge cases: an unsaved parameterized query, a name with no letter", () => {
  const session = feature(test, "features/queries/query-editor-edge-cases.feature", import.meta.url);
  test("Run query... on an unsaved parameterized query opens its parameter dialog", {tag: ["@edge-cases", "@realizes:views.queries", "@known-failure:43", "@realizes:core#641"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And no query named \"BDDQUnsaved{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDQUnsaved{time}")));
    await session.step(22, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(23, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await knownFailure(async () => {
      await session.step(30, "When user picks \"New Query...\" from the context menu of Databases---Postgres---Datagrok tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(31, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(32, "When user enters \"BDDQUnsaved{time}\" into Name input", () => enterInto(page, session.text("BDDQUnsaved{time}"), el("Name input")));
      await session.step(33, "And user replaces the code of code editor with '--input: string typeName = \"Project\"'", () => replaceCode(page, el("code editor"), "--input: string typeName = \"Project\""));
      await session.step(34, "And user appends \"select name from entity_types where name = @typeName\" to code editor", () => appendToEditor(page, "select name from entity_types where name = @typeName", el("code editor")));
      await session.step(35, "And user clicks on play icon", () => clickOn(page, el("play icon")));
      await session.step(36, "Then \"BDDQUnsaved{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDQUnsaved{time}\" dialog")), "visible"));
      await session.step(37, "When user clicks on OK button in \"BDDQUnsaved{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BDDQUnsaved{time}\" dialog"))));
      await session.step(38, "Then the \"BDDQUnsaved{time}\" dialog should close", () => dialogCloses(page, session.text("BDDQUnsaved{time}")));
      await session.step(39, "And grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(40, "And the \"rows\" reading of grid should be 1", () => readingIs(page, "rows", el("grid"), 1));
      await session.step(41, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(42, "When user clicks on \"Run query...\" action in toolbox", () => clickOn(page, el("\"Run query...\" action in toolbox")));
      await session.step(43, "Then \"BDDQUnsaved{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDQUnsaved{time}\" dialog")), "visible"));
      await session.step(44, "When user clicks on CANCEL button in \"BDDQUnsaved{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDDQUnsaved{time}\" dialog"))));
      await session.step(45, "Then the \"BDDQUnsaved{time}\" dialog should close", () => dialogCloses(page, session.text("BDDQUnsaved{time}")));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    }, 43);
  });
  test("Rename... refuses a name with no letter, and CANCEL keeps the name", {tag: ["@edge-cases", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And no query named \"BDDQUnsaved{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDQUnsaved{time}")));
    await session.step(22, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(23, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(49, "Given a query \"BDDQRenameEdge{time}\" on \"System:Datagrok\" reads \"select name from entity_types\"", () => queryOnServer(page, session.text("BDDQRenameEdge{time}"), "System:Datagrok", "select name from entity_types"));
    await session.step(50, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(51, "When user picks \"Rename...\" from the context menu of Databases---Postgres---Datagrok---BDDQRenameEdge{time} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Databases---Postgres---Datagrok---BDDQRenameEdge{time} tree node inside browse tree"))));
    await session.step(52, "Then \"Rename dataquery\" dialog should be visible", () => shouldBe(page, el("\"Rename dataquery\" dialog"), "visible"));
    await session.step(53, "When user enters \"123\" into Name input in \"Rename dataquery\" dialog", () => enterInto(page, "123", el("Name input in \"Rename dataquery\" dialog")));
    await session.step(54, "Then OK button in \"Rename dataquery\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Rename dataquery\" dialog"), "disabled"));
    await session.step(55, "When user enters \"   \" into Name input in \"Rename dataquery\" dialog", () => enterInto(page, "   ", el("Name input in \"Rename dataquery\" dialog")));
    await session.step(56, "Then OK button in \"Rename dataquery\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Rename dataquery\" dialog"), "disabled"));
    await session.step(57, "When user enters \"BDDQRenamedEdge{time}\" into Name input in \"Rename dataquery\" dialog", () => enterInto(page, session.text("BDDQRenamedEdge{time}"), el("Name input in \"Rename dataquery\" dialog")));
    await session.step(58, "Then OK button in \"Rename dataquery\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Rename dataquery\" dialog"), "enabled"));
    await session.step(59, "When user clicks on CANCEL button in \"Rename dataquery\" dialog", () => clickOn(page, el("CANCEL button in \"Rename dataquery\" dialog")));
    await session.step(60, "Then the \"Rename dataquery\" dialog should close", () => dialogCloses(page, "Rename dataquery"));
    await session.step(61, "And 1 query named \"BDDQRenameEdge{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDDQRenameEdge{time}")));
    await session.step(62, "And 0 queries named \"BDDQRenamedEdge{time}\" should be on the server", () => queriesOnServer(page, 0, session.text("BDDQRenamedEdge{time}")));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
  });
});
