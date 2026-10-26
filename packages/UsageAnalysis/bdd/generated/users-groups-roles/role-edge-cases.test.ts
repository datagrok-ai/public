/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/role-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.roles]
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
import {check, clickOn, expand, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, galleryCountHigher, galleryCountLower, groupsOnServer, noGroupOnServer, rememberGalleryCount, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Role edge cases: a name a group already has, a name of spaces only, a cancelled grant", () => {
  const session = feature(test, "features/users-groups-roles/role-edge-cases.feature", import.meta.url);
  test("A new role under the name of a group is refused", {tag: ["@serial", "@roles", "@edge-cases", "@realizes:views.roles"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And no group named \"BDD-RE-Group-{time}, BDD-RE-Role-{time}\" is on the server", () => noGroupOnServer(page, session.text("BDD-RE-Group-{time}, BDD-RE-Role-{time}")));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(22, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(25, "When user clicks on \"Platform > Groups\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Groups\" tree node inside browse tree")));
    await session.step(26, "Then the \"Groups\" view should be current", () => viewIsCurrent(page, "Groups"));
    await session.step(27, "When user clicks on \"New Group...\" button", () => clickOn(page, el("\"New Group...\" button")));
    await session.step(28, "And user types \"BDD-RE-Group-{time}\" into Name input in \"Create New Group\" dialog", () => typeInto(page, session.text("BDD-RE-Group-{time}"), el("Name input in \"Create New Group\" dialog")));
    await session.step(29, "And user clicks on OK button in \"Create New Group\" dialog", () => clickOn(page, el("OK button in \"Create New Group\" dialog")));
    await session.step(30, "Then the \"Create New Group\" dialog should close", () => dialogCloses(page, "Create New Group"));
    await session.step(31, "And 1 group named \"BDD-RE-Group-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-RE-Group-{time}")));
    await session.step(32, "When user clicks on \"Platform > Roles\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Roles\" tree node inside browse tree")));
    await session.step(33, "Then the \"Roles\" view should be current", () => viewIsCurrent(page, "Roles"));
    await session.step(34, "When user clicks on \"New Role...\" button", () => clickOn(page, el("\"New Role...\" button")));
    await session.step(35, "And user types \"BDD-RE-Group-{time}\" into Name input in \"Create New Role\" dialog", () => typeInto(page, session.text("BDD-RE-Group-{time}"), el("Name input in \"Create New Role\" dialog")));
    await session.step(36, "And user clicks on OK button in \"Create New Role\" dialog", () => clickOn(page, el("OK button in \"Create New Role\" dialog")));
    await session.step(37, "Then an error balloon containing \"A group or role with the same name already exists\" should have been shown", () => errorBalloonText(page, "A group or role with the same name already exists"));
    await session.step(38, "And \"Create New Role\" dialog should be visible", () => shouldBe(page, el("\"Create New Role\" dialog"), "visible"));
    await session.step(39, "And 1 group named \"BDD-RE-Group-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-RE-Group-{time}")));
    await session.step(40, "When user clicks on CANCEL button in \"Create New Role\" dialog", () => clickOn(page, el("CANCEL button in \"Create New Role\" dialog")));
    await session.step(41, "Then the \"Create New Role\" dialog should close", () => dialogCloses(page, "Create New Role"));
  });
  test("The New Role dialog refuses a name of spaces only", {tag: ["@serial", "@roles", "@edge-cases", "@realizes:views.roles"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And no group named \"BDD-RE-Group-{time}, BDD-RE-Role-{time}\" is on the server", () => noGroupOnServer(page, session.text("BDD-RE-Group-{time}, BDD-RE-Role-{time}")));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(22, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(44, "When user clicks on \"Platform > Roles\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Roles\" tree node inside browse tree")));
    await session.step(45, "Then the \"Roles\" view should be current", () => viewIsCurrent(page, "Roles"));
    await session.step(46, "When user clicks on \"New Role...\" button", () => clickOn(page, el("\"New Role...\" button")));
    await session.step(47, "And user types \"   \" into Name input in \"Create New Role\" dialog", () => typeInto(page, "   ", el("Name input in \"Create New Role\" dialog")));
    await session.step(48, "Then Name input in \"Create New Role\" dialog should be invalid", () => shouldBe(page, el("Name input in \"Create New Role\" dialog"), "invalid"));
    await session.step(49, "And OK button in \"Create New Role\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create New Role\" dialog"), "disabled"));
    await session.step(50, "When user clicks on CANCEL button in \"Create New Role\" dialog", () => clickOn(page, el("CANCEL button in \"Create New Role\" dialog")));
    await session.step(51, "Then the \"Create New Role\" dialog should close", () => dialogCloses(page, "Create New Role"));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
    await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("CANCEL in Global Permissions grants nothing", {tag: ["@serial", "@roles", "@edge-cases", "@realizes:views.roles"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And no group named \"BDD-RE-Group-{time}, BDD-RE-Role-{time}\" is on the server", () => noGroupOnServer(page, session.text("BDD-RE-Group-{time}, BDD-RE-Role-{time}")));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(22, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(56, "When user clicks on \"Platform > Roles\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Roles\" tree node inside browse tree")));
    await session.step(57, "Then the \"Roles\" view should be current", () => viewIsCurrent(page, "Roles"));
    await session.step(58, "When user remembers the gallery counter", () => rememberGalleryCount(page));
    await session.step(59, "And user clicks on \"New Role...\" button", () => clickOn(page, el("\"New Role...\" button")));
    await session.step(60, "And user types \"BDD-RE-Role-{time}\" into Name input in \"Create New Role\" dialog", () => typeInto(page, session.text("BDD-RE-Role-{time}"), el("Name input in \"Create New Role\" dialog")));
    await session.step(61, "And user clicks on OK button in \"Create New Role\" dialog", () => clickOn(page, el("OK button in \"Create New Role\" dialog")));
    await session.step(62, "Then the \"Create New Role\" dialog should close", () => dialogCloses(page, "Create New Role"));
    await session.step(63, "And 1 role named \"BDD-RE-Role-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-RE-Role-{time}")));
    await session.step(64, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
    await session.step(65, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
    await session.step(66, "When user remembers the gallery counter", () => rememberGalleryCount(page));
    await session.step(67, "And user types \"BDD-RE-Role-{time}\" into gallery search", () => typeInto(page, session.text("BDD-RE-Role-{time}"), el("gallery search")));
    await session.step(68, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
    await session.step(69, "When user clicks on \"BDD-RE-Role-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-RE-Role-{time}\" link in gallery"))));
    await session.step(70, "Then the context panel should show \"BDD-RE-Role-{time}\"", () => contextPanelShows(page, session.text("BDD-RE-Role-{time}")));
    await session.step(71, "When user expands \"Global Permissions\" section in context panel", () => expand(page, el("\"Global Permissions\" section in context panel")));
    await session.step(72, "And user clicks on MANAGE button in \"Global Permissions\" section in context panel", () => clickOn(page, el("MANAGE button in \"Global Permissions\" section in context panel")));
    await session.step(73, "Then \"BDD-RE-Role-{time}: Global Permissions\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-RE-Role-{time}: Global Permissions\" dialog")), "visible"));
    await session.step(74, "When user expands \"Browse\" tree node in \"BDD-RE-Role-{time}: Global Permissions\" dialog", () => expand(page, el(session.text("\"Browse\" tree node in \"BDD-RE-Role-{time}: Global Permissions\" dialog"))));
    await session.step(75, "And user checks \"Browse > Browse Apps\" tree node in \"BDD-RE-Role-{time}: Global Permissions\" dialog", () => check(page, el(session.text("\"Browse > Browse Apps\" tree node in \"BDD-RE-Role-{time}: Global Permissions\" dialog"))));
    await session.step(76, "Then \"Browse > Browse Apps\" tree node in \"BDD-RE-Role-{time}: Global Permissions\" dialog should be checked", () => shouldBe(page, el(session.text("\"Browse > Browse Apps\" tree node in \"BDD-RE-Role-{time}: Global Permissions\" dialog")), "checked"));
    await session.step(77, "When user clicks on CANCEL button in \"BDD-RE-Role-{time}: Global Permissions\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-RE-Role-{time}: Global Permissions\" dialog"))));
    await session.step(78, "Then the \"BDD-RE-Role-{time}: Global Permissions\" dialog should close", () => dialogCloses(page, session.text("BDD-RE-Role-{time}: Global Permissions")));
    await session.step(79, "When user clicks on MANAGE button in \"Global Permissions\" section in context panel", () => clickOn(page, el("MANAGE button in \"Global Permissions\" section in context panel")));
    await session.step(80, "And user expands \"Browse\" tree node in \"BDD-RE-Role-{time}: Global Permissions\" dialog", () => expand(page, el(session.text("\"Browse\" tree node in \"BDD-RE-Role-{time}: Global Permissions\" dialog"))));
    await session.step(81, "Then \"Browse > Browse Apps\" tree node in \"BDD-RE-Role-{time}: Global Permissions\" dialog should be unchecked", () => shouldBe(page, el(session.text("\"Browse > Browse Apps\" tree node in \"BDD-RE-Role-{time}: Global Permissions\" dialog")), "unchecked"));
    await session.step(82, "When user clicks on CANCEL button in \"BDD-RE-Role-{time}: Global Permissions\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-RE-Role-{time}: Global Permissions\" dialog"))));
    await session.step(83, "Then the \"BDD-RE-Role-{time}: Global Permissions\" dialog should close", () => dialogCloses(page, session.text("BDD-RE-Role-{time}: Global Permissions")));
    await session.step(84, "And \"Global Permissions\" section in context panel should contain text \"No global permissions\"", () => shouldContainText(page, el("\"Global Permissions\" section in context panel"), "No global permissions"));
    await session.step(85, "And no errors should have been logged", () => noErrors(page));
    await session.step(86, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
