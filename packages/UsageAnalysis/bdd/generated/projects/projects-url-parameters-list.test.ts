/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-url-parameters-list.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, views.functions, core#567]
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
import {clickOn, enterInto, hoverOver, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {dblclickRunsQuery} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, contextPanelOpen, dialogCloses, noProjectOnServer, queryTextOnServer, simpleModeOff, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pointerAway} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("A dashboard on a query with a list parameter: the Source link after the save", () => {
  const session = feature(test, "features/projects/projects-url-parameters-list.feature", import.meta.url);
  test("After the save, the Source link still carries the scalar parameter", {tag: ["@realizes:views.projects", "@realizes:views.functions", "@known-failure:54", "@realizes:core#567"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And no project named \"BDDUrlListProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDUrlListProj{time}")));
    await session.step(22, "And a query \"BDDUrlListQuery{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDUrlListQuery{time}"), "System:Datagrok", "--input: string typeName = \"Project\"\n--input: list<string> otherNames\nselect name from entity_types where name = @typeName or name = any(@otherNames)"));
    await knownFailure(async () => {
      await session.step(32, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(33, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(34, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(35, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(36, "When user double-clicks on Databases---Postgres---Datagrok---BDDUrlListQuery{time} tree node inside browse tree and the query it runs completes", () => dblclickRunsQuery(page, el(session.text("Databases---Postgres---Datagrok---BDDUrlListQuery{time} tree node inside browse tree"))));
      await session.step(37, "Then the \"BDDUrlListQuery{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlListQuery{time}")));
      await session.step(38, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(39, "When user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(40, "Then tooltip should contain text \"/func/\"", () => shouldContainText(page, el("tooltip"), "/func/"));
      await session.step(41, "And tooltip should contain text \"typeName=%22Project%22\"", () => shouldContainText(page, el("tooltip"), "typeName=%22Project%22"));
      await session.step(42, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(43, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(44, "When user enters \"BDDUrlListProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDUrlListProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(45, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(46, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(47, "And \"Share BDDUrlListProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDUrlListProj{time}\" dialog")), "visible"));
      await session.step(48, "When user clicks on CANCEL button in \"Share BDDUrlListProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDUrlListProj{time}\" dialog"))));
      await session.step(49, "Then the \"Share BDDUrlListProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDUrlListProj{time}")));
      await session.step(50, "When user moves the pointer away from Source pane in toolbox", () => pointerAway(page, el("Source pane in toolbox")));
      await session.step(51, "And user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(52, "Then tooltip should contain text \"/p/\"", () => shouldContainText(page, el("tooltip"), "/p/"));
      await session.step(53, "And tooltip should contain text \".BDDUrlListProj{time}\"", () => shouldContainText(page, el("tooltip"), session.text(".BDDUrlListProj{time}")));
      await session.step(54, "And tooltip should contain text \"typeName=Project\"", () => shouldContainText(page, el("tooltip"), "typeName=Project"));
    }, 54);
  });
});
