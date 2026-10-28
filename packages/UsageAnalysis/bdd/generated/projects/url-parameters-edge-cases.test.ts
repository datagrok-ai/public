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
import {clickOn, enterInto, isExpanded, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {dblclickRunsQuery} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, contextPanelOpen, dialogCloses, noProjectOnServer, queryTextOnServer, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorOrWarningBalloonMatching, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("URL parameters edge cases: a value of the wrong type", () => {
  const session = feature(test, "features/projects/url-parameters-edge-cases.feature", import.meta.url);
  test("A value of the wrong type in the address is refused with a balloon naming the parameter", {tag: ["@realizes:views.projects", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And no project named \"BDDUrlEdgeP{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDUrlEdgeP{time}")));
    await session.step(20, "And a query \"BDDUrlEdgeQ{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDUrlEdgeQ{time}"), "System:Datagrok", "--input: int lim = 2\n--input: string typeName = \"Project\"\nselect name from entity_types where name <> @typeName order by name limit @lim"));
    await session.step(26, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
    await session.step(27, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(28, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(29, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(30, "When user double-clicks on Databases---Postgres---Datagrok---BDDUrlEdgeQ{time} tree node inside browse tree and the query it runs completes", () => dblclickRunsQuery(page, el(session.text("Databases---Postgres---Datagrok---BDDUrlEdgeQ{time} tree node inside browse tree"))));
    await session.step(31, "Then the \"BDDUrlEdgeQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlEdgeQ{time}")));
    await session.step(32, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(33, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(34, "When user enters \"BDDUrlEdgeP{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDUrlEdgeP{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(35, "And user clicks on \"URL Parameters\" button in \"BDDUrlEdgeQ{time}\" project table in \"Save project\" dialog", () => clickOn(page, el(session.text("\"URL Parameters\" button in \"BDDUrlEdgeQ{time}\" project table in \"Save project\" dialog"))));
    await session.step(36, "Then first \"URL alias\" text input in \"Save project\" dialog should have value \"lim\"", () => shouldHaveValue(page, el("first \"URL alias\" text input in \"Save project\" dialog"), "lim"));
    await session.step(40, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(41, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(42, "And \"Share BDDUrlEdgeP{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDUrlEdgeP{time}\" dialog")), "visible"));
    await session.step(43, "When user clicks on CANCEL button in \"Share BDDUrlEdgeP{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDUrlEdgeP{time}\" dialog"))));
    await session.step(44, "Then the \"Share BDDUrlEdgeP{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDUrlEdgeP{time}")));
    await session.step(45, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
    await session.step(46, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    await session.step(47, "When user opens the address \"/p/Admin.BDDUrlEdgeP{time}?lim=abc\"", () => openAddress(page, session.text("/p/Admin.BDDUrlEdgeP{time}?lim=abc")));
    await session.step(48, "Then an error or warning balloon matching \"lim\" should have been shown", () => errorOrWarningBalloonMatching(page, "lim"));
  });
});
