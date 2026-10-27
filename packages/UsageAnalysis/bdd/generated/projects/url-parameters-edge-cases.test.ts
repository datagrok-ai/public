/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/url-parameters-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects]
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
import {clearField, clickOn, enterInto, isExpanded, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {dblclickRunsQuery} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, contextPanelOpen, dialogCloses, noProjectOnServer, projectsOnServer, queryTextOnServer, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorOrWarningBalloonMatching, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("URL parameters edge cases: a value of the wrong type, an emptied alias, one alias for two parameters", () => {
  const session = feature(test, "features/projects/url-parameters-edge-cases.feature", import.meta.url);
  test("A value of the wrong type in the address is refused with a balloon naming the parameter", {tag: ["@realizes:views.projects", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(23, "And simple mode is off", () => simpleModeOff(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And no project named \"BDDUrlEdgeP{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDUrlEdgeP{time}")));
    await session.step(26, "And a query \"BDDUrlEdgeQ{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDUrlEdgeQ{time}"), "System:Datagrok", "--input: int lim = 2\n--input: string typeName = \"Project\"\nselect name from entity_types where name <> @typeName order by name limit @lim"));
    await session.step(32, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
    await session.step(33, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(34, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(35, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(36, "When user double-clicks on Databases---Postgres---Datagrok---BDDUrlEdgeQ{time} tree node inside browse tree and the query it runs completes", () => dblclickRunsQuery(page, el(session.text("Databases---Postgres---Datagrok---BDDUrlEdgeQ{time} tree node inside browse tree"))));
    await session.step(37, "Then the \"BDDUrlEdgeQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlEdgeQ{time}")));
    await session.step(38, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(39, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(40, "When user enters \"BDDUrlEdgeP{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDUrlEdgeP{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(41, "And user clicks on \"URL Parameters\" button in \"BDDUrlEdgeQ{time}\" project table in \"Save project\" dialog", () => clickOn(page, el(session.text("\"URL Parameters\" button in \"BDDUrlEdgeQ{time}\" project table in \"Save project\" dialog"))));
    await session.step(42, "Then first \"URL alias\" text input in \"Save project\" dialog should have value \"lim\"", () => shouldHaveValue(page, el("first \"URL alias\" text input in \"Save project\" dialog"), "lim"));
    await session.step(46, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(47, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(48, "And \"Share BDDUrlEdgeP{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDUrlEdgeP{time}\" dialog")), "visible"));
    await session.step(49, "When user clicks on CANCEL button in \"Share BDDUrlEdgeP{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDUrlEdgeP{time}\" dialog"))));
    await session.step(50, "Then the \"Share BDDUrlEdgeP{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDUrlEdgeP{time}")));
    await session.step(51, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
    await session.step(52, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    await session.step(53, "When user opens the address \"/p/Admin.BDDUrlEdgeP{time}?lim=abc\"", () => openAddress(page, session.text("/p/Admin.BDDUrlEdgeP{time}?lim=abc")));
    await session.step(54, "Then an error or warning balloon matching \"lim\" should have been shown", () => errorOrWarningBalloonMatching(page, "lim"));
  });
  test("An emptied alias falls back to the parameter's own name", {tag: ["@realizes:views.projects", "@edge-cases", "@known-failure:59", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(23, "And simple mode is off", () => simpleModeOff(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And no project named \"BDDUrlEdgeP{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDUrlEdgeP{time}")));
    await session.step(26, "And a query \"BDDUrlEdgeQ{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDUrlEdgeQ{time}"), "System:Datagrok", "--input: int lim = 2\n--input: string typeName = \"Project\"\nselect name from entity_types where name <> @typeName order by name limit @lim"));
    await session.step(32, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
    await session.step(33, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(34, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(35, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(36, "When user double-clicks on Databases---Postgres---Datagrok---BDDUrlEdgeQ{time} tree node inside browse tree and the query it runs completes", () => dblclickRunsQuery(page, el(session.text("Databases---Postgres---Datagrok---BDDUrlEdgeQ{time} tree node inside browse tree"))));
    await session.step(37, "Then the \"BDDUrlEdgeQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlEdgeQ{time}")));
    await session.step(38, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(39, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(40, "When user enters \"BDDUrlEdgeP{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDUrlEdgeP{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(41, "And user clicks on \"URL Parameters\" button in \"BDDUrlEdgeQ{time}\" project table in \"Save project\" dialog", () => clickOn(page, el(session.text("\"URL Parameters\" button in \"BDDUrlEdgeQ{time}\" project table in \"Save project\" dialog"))));
    await session.step(42, "Then first \"URL alias\" text input in \"Save project\" dialog should have value \"lim\"", () => shouldHaveValue(page, el("first \"URL alias\" text input in \"Save project\" dialog"), "lim"));
    await knownFailure(async () => {
      await session.step(58, "When user clears first \"URL alias\" text input in \"Save project\" dialog", () => clearField(page, el("first \"URL alias\" text input in \"Save project\" dialog")));
      await session.step(59, "Then no errors should have been logged", () => noErrors(page));
      await session.step(60, "And \"Save project\" dialog should contain text \"?lim=2\"", () => shouldContainText(page, el("\"Save project\" dialog"), "?lim=2"));
      await session.step(61, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(62, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(63, "And 0 projects named \"BDDUrlEdgeP{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDUrlEdgeP{time}")));
    }, 59);
  });
  test("Two parameters given one alias are refused when the dashboard is saved", {tag: ["@realizes:views.projects", "@edge-cases", "@known-failure:71", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(23, "And simple mode is off", () => simpleModeOff(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And no project named \"BDDUrlEdgeP{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDUrlEdgeP{time}")));
    await session.step(26, "And a query \"BDDUrlEdgeQ{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDUrlEdgeQ{time}"), "System:Datagrok", "--input: int lim = 2\n--input: string typeName = \"Project\"\nselect name from entity_types where name <> @typeName order by name limit @lim"));
    await session.step(32, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
    await session.step(33, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(34, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(35, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(36, "When user double-clicks on Databases---Postgres---Datagrok---BDDUrlEdgeQ{time} tree node inside browse tree and the query it runs completes", () => dblclickRunsQuery(page, el(session.text("Databases---Postgres---Datagrok---BDDUrlEdgeQ{time} tree node inside browse tree"))));
    await session.step(37, "Then the \"BDDUrlEdgeQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlEdgeQ{time}")));
    await session.step(38, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(39, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(40, "When user enters \"BDDUrlEdgeP{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDUrlEdgeP{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(41, "And user clicks on \"URL Parameters\" button in \"BDDUrlEdgeQ{time}\" project table in \"Save project\" dialog", () => clickOn(page, el(session.text("\"URL Parameters\" button in \"BDDUrlEdgeQ{time}\" project table in \"Save project\" dialog"))));
    await session.step(42, "Then first \"URL alias\" text input in \"Save project\" dialog should have value \"lim\"", () => shouldHaveValue(page, el("first \"URL alias\" text input in \"Save project\" dialog"), "lim"));
    await knownFailure(async () => {
      await session.step(69, "When user enters \"typeName\" into first \"URL alias\" text input in \"Save project\" dialog", () => enterInto(page, "typeName", el("first \"URL alias\" text input in \"Save project\" dialog")));
      await session.step(70, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(71, "Then an error or warning balloon matching \"alias\" should have been shown", () => errorOrWarningBalloonMatching(page, "alias"));
      await session.step(72, "And \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(73, "And 0 projects named \"BDDUrlEdgeP{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDUrlEdgeP{time}")));
      await session.step(74, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(75, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    }, 71);
  });
});
