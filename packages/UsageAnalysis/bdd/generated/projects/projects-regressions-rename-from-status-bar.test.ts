/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-rename-from-status-bar.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, core#626]
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
import {clickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {dblclickRunsQuery} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, queryTextOnServer, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a query result renamed from the status bar", () => {
  const session = feature(test, "features/projects/projects-regressions-rename-from-status-bar.feature", import.meta.url);
  test("The table renamed through the status bar takes the new name with its view", {tag: ["@realizes:views.projects", "@known-failure:38", "@realizes:core#626"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And a query \"BDDRenameQ{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDRenameQ{time}"), "System:Datagrok", "select name from entity_types"));
    await knownFailure(async () => {
      await session.step(24, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
      await session.step(25, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(26, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(27, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(28, "When user double-clicks on Databases---Postgres---Datagrok---BDDRenameQ{time} tree node inside browse tree and the query it runs completes", () => dblclickRunsQuery(page, el(session.text("Databases---Postgres---Datagrok---BDDRenameQ{time} tree node inside browse tree"))));
      await session.step(29, "Then the \"BDDRenameQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDRenameQ{time}")));
      await session.step(30, "When user clicks on \"BDDRenameQ{time}\" text in status bar", () => clickOn(page, el(session.text("\"BDDRenameQ{time}\" text in status bar"))));
      await session.step(31, "Then the context panel should show \"BDDRenameQ{time}\"", () => contextPanelShows(page, session.text("BDDRenameQ{time}")));
      await session.step(32, "Given Actions pane in context panel is expanded", () => isExpanded(page, el("Actions pane in context panel")));
      await session.step(33, "When user clicks on \"Rename...\" label in Actions pane in context panel", () => clickOn(page, el("\"Rename...\" label in Actions pane in context panel")));
      await session.step(34, "Then \"Rename table\" dialog should be visible", () => shouldBe(page, el("\"Rename table\" dialog"), "visible"));
      await session.step(35, "When user enters \"BDDRenamedT{time}\" into \"New name:\" text input in \"Rename table\" dialog", () => enterInto(page, session.text("BDDRenamedT{time}"), el("\"New name:\" text input in \"Rename table\" dialog")));
      await session.step(36, "And user clicks on OK button in \"Rename table\" dialog", () => clickOn(page, el("OK button in \"Rename table\" dialog")));
      await session.step(37, "Then the \"Rename table\" dialog should close", () => dialogCloses(page, "Rename table"));
      await session.step(38, "And the open table views should be exactly \"BDDRenamedT{time}\"", () => openTableViewsExactly(page, session.text("BDDRenamedT{time}")));
      await session.step(39, "And the \"BDDRenamedT{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDRenamedT{time}")));
    }, 38);
  });
});
