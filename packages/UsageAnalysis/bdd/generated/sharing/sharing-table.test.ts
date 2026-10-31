/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sharing/sharing-table.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog]
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
import {clickOn, expand, isExpanded, pressKey, shouldBe, shouldContainText, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, contextPanelOpen, contextPanelShows, dialogCloses, noProjectOnServer, openDatasetRowsAs, openSaveDialog, pickSharingUser, projectsOnServer, refreshBrowse, saveFromDialog, sharingUserLevel} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A table saved in a project, shared from its Sharing pane", () => {
  const session = feature(test, "features/sharing/sharing-table.feature", import.meta.url);
  test("Share... on the table shares it with another user", {tag: ["@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And no project named \"BDD-ShareTable-{time}\" is on the server", () => noProjectOnServer(page, session.text("BDD-ShareTable-{time}")));
    await session.step(24, "And user opens demog-1000 dataset keeping the first 20 rows as \"BDD-ShareTable-{time}\"", () => openDatasetRowsAs(page, ds("demog-1000"), 20, session.text("BDD-ShareTable-{time}")));
    await session.step(25, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
    await session.step(26, "And user types \"BDD-ShareTable-{time}\" into text input in \"Save project\" dialog", () => typeInto(page, session.text("BDD-ShareTable-{time}"), el("text input in \"Save project\" dialog")));
    await session.step(27, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
    await session.step(28, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(29, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(30, "Then 1 project named \"BDD-ShareTable-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDD-ShareTable-{time}")));
    await session.step(31, "When user closes all views", () => closeAllViews(page));
    await session.step(32, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(34, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
    await session.step(35, "And \"My stuff > My dashboards\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > My dashboards\" tree node inside browse tree")));
    await session.step(36, "And \"My stuff > My dashboards > BDD-ShareTable-{time}\" tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("\"My stuff > My dashboards > BDD-ShareTable-{time}\" tree node inside browse tree"))));
    await session.step(39, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(40, "When user clicks on \"My stuff > My dashboards > BDD-ShareTable-{time} > BDD-ShareTable-{time}\" tree node inside browse tree", () => clickOn(page, el(session.text("\"My stuff > My dashboards > BDD-ShareTable-{time} > BDD-ShareTable-{time}\" tree node inside browse tree"))));
    await session.step(41, "Then the context panel should show \"BDD-ShareTable-{time}\"", () => contextPanelShows(page, session.text("BDD-ShareTable-{time}")));
    await session.step(42, "When user expands Sharing section in context panel", () => expand(page, el("Sharing section in context panel")));
    await session.step(43, "Then Sharing section in context panel should contain text \"You are the owner\"", () => shouldContainText(page, el("Sharing section in context panel"), "You are the owner"));
    await session.step(44, "When user clicks on \"Share...\" button in Sharing section in context panel", () => clickOn(page, el("\"Share...\" button in Sharing section in context panel")));
    await session.step(45, "Then \"Share BDD-ShareTable-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD-ShareTable-{time}\" dialog")), "visible"));
    await session.step(46, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
    await session.step(47, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDD-ShareTable-{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDD-ShareTable-{time}\" dialog"))));
    await session.step(48, "And user unchecks \"Send notifications\" input in \"Share BDD-ShareTable-{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDD-ShareTable-{time}\" dialog"))));
    await session.step(49, "And user clicks on OK button in \"Share BDD-ShareTable-{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDD-ShareTable-{time}\" dialog"))));
    await session.step(50, "Then the \"Share BDD-ShareTable-{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDD-ShareTable-{time}")));
    await session.step(51, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
    await session.step(52, "When user clicks on \"Share...\" button in Sharing section in context panel", () => clickOn(page, el("\"Share...\" button in Sharing section in context panel")));
    await session.step(53, "Then \"Share BDD-ShareTable-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD-ShareTable-{time}\" dialog")), "visible"));
    await session.step(54, "And the access level of the sharing user in \"Share BDD-ShareTable-{time}\" dialog should be \"View and use\"", () => sharingUserLevel(page, el(session.text("\"Share BDD-ShareTable-{time}\" dialog")), "View and use"));
    await session.step(55, "When user clicks on CANCEL button in \"Share BDD-ShareTable-{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDD-ShareTable-{time}\" dialog"))));
    await session.step(56, "Then the \"Share BDD-ShareTable-{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDD-ShareTable-{time}")));
    await session.step(57, "And no errors should have been logged", () => noErrors(page));
    await session.step(58, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
