/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-query.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, core#463, GROK-20174]
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
import {appendToEditor, clickOn, collapse, doubleClickOn, enterInto, isExpanded, replaceCode, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, closeCurrentView, currentViewType, dialogCloses, noProjectOnServer, queryText, queryTextNot, queryTextOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: projects of a parameterized query", () => {
  const session = feature(test, "features/projects/projects-regressions-query.feature", import.meta.url);
  test("Projects regressions: projects of a parameterized query", {tag: ["@journey", "@realizes:views.projects", "@realizes:core#463", "@known-failure:76", "@realizes:GROK-20174"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("A project of a query is saved, and the query loses one of its parameters in its editor", async () => {
      await session.step(36, "Given no project named \"BDDRegParamLess{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegParamLess{time}")));
      await session.step(37, "And a query \"BDDRegParamsLess{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDRegParamsLess{time}"), "System:Datagrok", "--input: string typeName = \"Project\"\n--input: string pattern = \"%\"\nselect name from entity_types where name = @typeName and name like @pattern"));
      await session.step(43, "When user closes all views", () => closeAllViews(page));
      await session.step(44, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(45, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(46, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(47, "Given Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(48, "When user double-clicks Databases---Postgres---Datagrok---BDDRegParamsLess{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("Databases---Postgres---Datagrok---BDDRegParamsLess{time} tree node inside browse tree"))));
      await session.step(49, "Then the \"BDDRegParamsLess{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDRegParamsLess{time}")));
      await session.step(50, "And the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
      await session.step(51, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(52, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(53, "When user enters \"BDDRegParamLess{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegParamLess{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(54, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(55, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(56, "And \"Share BDDRegParamLess{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegParamLess{time}\" dialog")), "visible"));
      await session.step(57, "When user clicks on CANCEL button in \"Share BDDRegParamLess{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegParamLess{time}\" dialog"))));
      await session.step(58, "Then the \"Share BDDRegParamLess{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegParamLess{time}")));
      await session.step(59, "When user closes all views", () => closeAllViews(page));
      await session.step(60, "When user collapses Databases---Postgres---Datagrok---BDDRegParamsLess{time} tree node inside browse tree", () => collapse(page, el(session.text("Databases---Postgres---Datagrok---BDDRegParamsLess{time} tree node inside browse tree"))));
      await session.step(61, "And user picks \"Edit...\" from the context menu of Databases---Postgres---Datagrok---BDDRegParamsLess{time} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---Datagrok---BDDRegParamsLess{time} tree node inside browse tree"))));
      await session.step(62, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(63, "When user replaces the code of code editor with '--input: string typeName = \"Project\"'", () => replaceCode(page, el("code editor"), "--input: string typeName = \"Project\""));
      await session.step(64, "And user appends \"select name from entity_types where name = @typeName\" to code editor", () => appendToEditor(page, "select name from entity_types where name = @typeName", el("code editor")));
      await session.step(65, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(66, "Then the query \"BDDRegParamsLess{time}\" on the server should not read \"pattern\"", () => queryTextNot(page, session.text("BDDRegParamsLess{time}"), "pattern"));
      await session.step(67, "When user closes the current view", () => closeCurrentView(page));
      await session.step(68, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(69, "And user enters \"BDDRegParamLess{time}\" into gallery search", () => enterInto(page, session.text("BDDRegParamLess{time}"), el("gallery search")));
      await session.step(70, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(71, "Then BDDRegParamLess{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDRegParamLess{time} gallery card")), "visible"));
    });
    await run.scenario("The project opens after its query lost a parameter (core#463)", async () => {
      await session.step(75, "When user double-clicks on BDDRegParamLess{time} gallery card", () => doubleClickOn(page, el(session.text("BDDRegParamLess{time} gallery card"))));
      await session.step(76, "Then the \"BDDRegParamsLess{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDRegParamsLess{time}")));
      await session.step(77, "And the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
      await session.step(78, "And \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
    }, {knownFailure: true, line: 76});
    await run.scenario("A project whose query gained a parameter names the parameters, not \"related entities were deleted\"", async () => {
      await session.step(82, "Given no project named \"BDDRegParamMore{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegParamMore{time}")));
      await session.step(83, "And a query \"BDDRegParamsMore{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDRegParamsMore{time}"), "System:Datagrok", "--input: string typeName = \"Project\"\nselect name from entity_types where name = @typeName"));
      await session.step(88, "When user closes all views", () => closeAllViews(page));
      await session.step(89, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(90, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(91, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(92, "Given Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(93, "When user double-clicks Databases---Postgres---Datagrok---BDDRegParamsMore{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("Databases---Postgres---Datagrok---BDDRegParamsMore{time} tree node inside browse tree"))));
      await session.step(94, "Then the \"BDDRegParamsMore{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDRegParamsMore{time}")));
      await session.step(95, "And the \"text of cell 1 of name\" reading of grid should be \"Project\"", () => readingReads(page, "text of cell 1 of name", el("grid"), "Project"));
      await session.step(96, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(97, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(98, "When user enters \"BDDRegParamMore{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegParamMore{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(99, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(100, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(101, "And \"Share BDDRegParamMore{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegParamMore{time}\" dialog")), "visible"));
      await session.step(102, "When user clicks on CANCEL button in \"Share BDDRegParamMore{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegParamMore{time}\" dialog"))));
      await session.step(103, "Then the \"Share BDDRegParamMore{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegParamMore{time}")));
      await session.step(104, "When user closes all views", () => closeAllViews(page));
      await session.step(105, "When user collapses Databases---Postgres---Datagrok---BDDRegParamsMore{time} tree node inside browse tree", () => collapse(page, el(session.text("Databases---Postgres---Datagrok---BDDRegParamsMore{time} tree node inside browse tree"))));
      await session.step(106, "And user picks \"Edit...\" from the context menu of Databases---Postgres---Datagrok---BDDRegParamsMore{time} tree node inside browse tree", () => pickFromContextMenu(page, "Edit...", el(session.text("Databases---Postgres---Datagrok---BDDRegParamsMore{time} tree node inside browse tree"))));
      await session.step(107, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
      await session.step(108, "When user replaces the code of code editor with '--input: string pattern = \"%\"'", () => replaceCode(page, el("code editor"), "--input: string pattern = \"%\""));
      await session.step(109, "And user appends '--input: string typeName = \"Project\"' to code editor", () => appendToEditor(page, "--input: string typeName = \"Project\"", el("code editor")));
      await session.step(110, "And user appends \"select name from entity_types where name = @typeName and name like @pattern\" to code editor", () => appendToEditor(page, "select name from entity_types where name = @typeName and name like @pattern", el("code editor")));
      await session.step(111, "And user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(112, "Then the query \"BDDRegParamsMore{time}\" on the server should read '--input: string pattern = \"%\"'", () => queryText(page, session.text("BDDRegParamsMore{time}"), "--input: string pattern = \"%\""));
      await session.step(113, "When user closes the current view", () => closeCurrentView(page));
      await session.step(114, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(115, "And user enters \"BDDRegParamMore{time}\" into gallery search", () => enterInto(page, session.text("BDDRegParamMore{time}"), el("gallery search")));
      await session.step(116, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(117, "And user double-clicks on BDDRegParamMore{time} gallery card", () => doubleClickOn(page, el(session.text("BDDRegParamMore{time} gallery card"))));
      await session.step(118, "Then \"Data loading error\" dialog should be visible", () => shouldBe(page, el("\"Data loading error\" dialog"), "visible"));
      await session.step(119, "And \"Data loading error\" dialog should contain text \"BDDRegParamsMore{time}\"", () => shouldContainText(page, el("\"Data loading error\" dialog"), session.text("BDDRegParamsMore{time}")));
      await session.step(120, "And \"Data loading error\" dialog should contain text \"required input parameters\"", () => shouldContainText(page, el("\"Data loading error\" dialog"), "required input parameters"));
      await session.step(121, "And \"Data loading error\" dialog should not contain text \"related entities\"", () => shouldNotContainText(page, el("\"Data loading error\" dialog"), "related entities"));
      await session.step(122, "When user clicks on \"CLOSE PROJECT\" button in \"Data loading error\" dialog", () => clickOn(page, el("\"CLOSE PROJECT\" button in \"Data loading error\" dialog")));
      await session.step(123, "Then the \"Data loading error\" dialog should close", () => dialogCloses(page, "Data loading error"));
    });
    run.finish();
  });
});
