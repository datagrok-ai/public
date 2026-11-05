/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-move-same-name.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, views.space, core#563]
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
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {dblclickRunsQuery} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, contextPanelOpen, dialogCloses, doubleClickOpensProject, noProjectOnServer, noSpaceOnServer, queryTextOnServer, savedWithDataSync, simpleModeOff, spaceOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a dashboard named like its query, moved into a space", () => {
  const session = feature(test, "features/projects/projects-regressions-move-same-name.feature", import.meta.url);
  test("The dashboard named like its query opens after the move into a space", {tag: ["@realizes:views.projects", "@realizes:views.space", "@known-failure:67", "@realizes:core#563"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "And no project named \"BDDSameName{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSameName{time}")));
    await session.step(21, "And no space named \"BDDSameNameSpace{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDSameNameSpace{time}")));
    await session.step(22, "And a space named \"BDDSameNameSpace{time}\" is on the server", () => spaceOnServer(page, session.text("BDDSameNameSpace{time}")));
    await session.step(23, "And a query \"BDDSameName{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDSameName{time}"), "System:Datagrok", "--input: string typeName = \"Project\"\nselect name from entity_types where name = @typeName"));
    await knownFailure(async () => {
      await session.step(32, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(33, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(34, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(35, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(36, "When user double-clicks on Databases---Postgres---Datagrok---BDDSameName{time} tree node inside browse tree and the query it runs completes", () => dblclickRunsQuery(page, el(session.text("Databases---Postgres---Datagrok---BDDSameName{time} tree node inside browse tree"))));
      await session.step(37, "Then the \"BDDSameName{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDSameName{time}")));
      await session.step(38, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(39, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(40, "And Data sync switch in \"BDDSameName{time}\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el(session.text("Data sync switch in \"BDDSameName{time}\" project table in \"Save project\" dialog")), "checked"));
      await session.step(41, "And Name text input in \"Save project\" dialog should have value \"BDDSameName{time}\"", () => shouldHaveValue(page, el("Name text input in \"Save project\" dialog"), session.text("BDDSameName{time}")));
      await session.step(42, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(43, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(44, "And \"Share BDDSameName{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDSameName{time}\" dialog")), "visible"));
      await session.step(45, "When user clicks on CANCEL button in \"Share BDDSameName{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDSameName{time}\" dialog"))));
      await session.step(46, "Then the \"Share BDDSameName{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDSameName{time}")));
      await session.step(47, "And the \"BDDSameName{time}\" table of the \"BDDSameName{time}\" project should be saved with data sync", () => savedWithDataSync(page, session.text("BDDSameName{time}"), session.text("BDDSameName{time}")));
      await session.step(48, "When user closes all views", () => closeAllViews(page));
      await session.step(49, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(50, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(51, "And user enters \"BDDSameName{time}\" into gallery search", () => enterInto(page, session.text("BDDSameName{time}"), el("gallery search")));
      await session.step(52, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(53, "And user picks \"Move to Space...\" from the context menu of BDDSameName{time} gallery card", () => pickFromContextMenu(page, "Move to Space...", el(session.text("BDDSameName{time} gallery card"))));
      await session.step(54, "Then \"Move to space\" dialog should be visible", () => shouldBe(page, el("\"Move to space\" dialog"), "visible"));
      await session.step(55, "When user selects \"BDDSameNameSpace{time}\" in Space input in \"Move to space\" dialog", () => selectIn(page, session.text("BDDSameNameSpace{time}"), el("Space input in \"Move to space\" dialog")));
      await session.step(56, "And user clicks on OK button in \"Move to space\" dialog", () => clickOn(page, el("OK button in \"Move to space\" dialog")));
      await session.step(57, "Then the \"Move to space\" dialog should close", () => dialogCloses(page, "Move to space"));
      await session.step(58, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(59, "When user double-clicks on BDDSameNameSpace{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDSameNameSpace{time} tree node inside browse tree"))));
      await session.step(60, "Then the \"BDDSameNameSpace{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDSameNameSpace{time}")));
      await session.step(61, "And \"BDDSameName{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDSameName{time}\" link in gallery")), "visible"));
      await session.step(62, "When user reloads the page", () => reloadPage(page));
      await session.step(63, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(64, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(65, "And user enters \"BDDSameName{time}\" into gallery search", () => enterInto(page, session.text("BDDSameName{time}"), el("gallery search")));
      await session.step(66, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(67, "And user double-clicks on BDDSameName{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDSameName{time} gallery card"))));
      await session.step(68, "Then the \"BDDSameName{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDSameName{time}")));
      await session.step(69, "And \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
      await session.step(70, "And the table should have 1 rows", () => rowCount(page, 1));
    }, 67);
  });
});
