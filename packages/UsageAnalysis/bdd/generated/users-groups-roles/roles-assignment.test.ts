/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/roles-assignment.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.roles]
--- */
import {test} from '@playwright/test';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, expand, shouldBe, shouldContainText, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {adminMemberOnServer, browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, galleryCountLower, newUserOnServer, noRoleOnServer, notMemberOnServer, plainMemberOnServer, rememberGalleryCount, rolesOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Who holds a role, and what it grants", () => {
  const session = feature(test, "features/users-groups-roles/roles-assignment.feature", import.meta.url);
  test("Who holds a role, and what it grants", {tag: ["@journey", "@roles", "@realizes:views.roles", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And a new user \"opavlenko{time}r\" with email \"opavlenko+{time}r@datagrok.ai\" is on the server", () => newUserOnServer(page, session.text("opavlenko{time}r"), session.text("opavlenko+{time}r@datagrok.ai")));
    await session.step(19, "And no role named \"BDD-RA-Role-{time}\" is on the server", () => noRoleOnServer(page, session.text("BDD-RA-Role-{time}")));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(22, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(23, "And user clicks on \"Platform > Roles\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Roles\" tree node inside browse tree")));
    await run.scenario("A role to assign", async () => {
      await session.step(26, "Then the \"Roles\" view should be current", () => viewIsCurrent(page, "Roles"));
      await session.step(27, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(28, "And user clicks on \"New Role...\" button", () => clickOn(page, el("\"New Role...\" button")));
      await session.step(29, "And user types \"BDD-RA-Role-{time}\" into Name input in \"Create New Role\" dialog", () => typeInto(page, session.text("BDD-RA-Role-{time}"), el("Name input in \"Create New Role\" dialog")));
      await session.step(30, "And user clicks on OK button in \"Create New Role\" dialog", () => clickOn(page, el("OK button in \"Create New Role\" dialog")));
      await session.step(31, "Then the \"Create New Role\" dialog should close", () => dialogCloses(page, "Create New Role"));
      await session.step(32, "And 1 role named \"BDD-RA-Role-{time}\" should be on the server", () => rolesOnServer(page, 1, session.text("BDD-RA-Role-{time}")));
      await session.step(33, "When user types \"BDD-RA-Role-{time}\" into gallery search", () => typeInto(page, session.text("BDD-RA-Role-{time}"), el("gallery search")));
      await session.step(34, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(35, "When user clicks on \"BDD-RA-Role-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-RA-Role-{time}\" link in gallery"))));
      await session.step(36, "Then the context panel should show \"BDD-RA-Role-{time}\"", () => contextPanelShows(page, session.text("BDD-RA-Role-{time}")));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
      await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("MANAGE assigns the role to a user (Roles-11)", async () => {
      await session.step(41, "When user clicks on MANAGE button in \"Assigned to\" section in context panel", () => clickOn(page, el("MANAGE button in \"Assigned to\" section in context panel")));
      await session.step(42, "Then \"BDD-RA-Role-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-RA-Role-{time} members\" dialog")), "visible"));
      await session.step(43, "When user types \"opavlenko{time}r\" into membership search", () => typeInto(page, session.text("opavlenko{time}r"), el("membership search")));
      await session.step(44, "And user clicks on add button of \"opavlenko{time}r\" membership candidate", () => clickOn(page, el(session.text("add button of \"opavlenko{time}r\" membership candidate"))));
      await session.step(45, "Then \"opavlenko{time}r\" membership row should be visible", () => shouldBe(page, el(session.text("\"opavlenko{time}r\" membership row")), "visible"));
      await session.step(46, "And checkbox label of \"opavlenko{time}r\" membership row should have text \"Can assign\"", () => shouldHaveText(page, el(session.text("checkbox label of \"opavlenko{time}r\" membership row")), "Can assign"));
      await session.step(47, "When user clicks on SAVE button in \"BDD-RA-Role-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-RA-Role-{time} members\" dialog"))));
      await session.step(48, "Then the \"BDD-RA-Role-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time} members")));
      await session.step(49, "And \"opavlenko{time}r\" should be a plain member of \"BDD-RA-Role-{time}\" on the server", () => plainMemberOnServer(page, session.text("opavlenko{time}r"), session.text("BDD-RA-Role-{time}")));
      await session.step(50, "And \"Assigned to\" section in context panel should contain text \"opavlenko{time}r\"", () => shouldContainText(page, el("\"Assigned to\" section in context panel"), session.text("opavlenko{time}r")));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
      await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Can assign is saved for the assignee (Roles-13)", async () => {
      await session.step(55, "When user clicks on MANAGE button in \"Assigned to\" section in context panel", () => clickOn(page, el("MANAGE button in \"Assigned to\" section in context panel")));
      await session.step(56, "Then checkbox of \"opavlenko{time}r\" membership row should be unchecked", () => shouldBe(page, el(session.text("checkbox of \"opavlenko{time}r\" membership row")), "unchecked"));
      await session.step(57, "When user checks checkbox of \"opavlenko{time}r\" membership row", () => check(page, el(session.text("checkbox of \"opavlenko{time}r\" membership row"))));
      await session.step(58, "And user clicks on SAVE button in \"BDD-RA-Role-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-RA-Role-{time} members\" dialog"))));
      await session.step(59, "Then the \"BDD-RA-Role-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time} members")));
      await session.step(60, "And \"opavlenko{time}r\" should be an admin member of \"BDD-RA-Role-{time}\" on the server", () => adminMemberOnServer(page, session.text("opavlenko{time}r"), session.text("BDD-RA-Role-{time}")));
      await session.step(61, "When user clicks on MANAGE button in \"Assigned to\" section in context panel", () => clickOn(page, el("MANAGE button in \"Assigned to\" section in context panel")));
      await session.step(62, "Then checkbox of \"opavlenko{time}r\" membership row should be checked", () => shouldBe(page, el(session.text("checkbox of \"opavlenko{time}r\" membership row")), "checked"));
      await session.step(63, "When user clicks on CANCEL button in \"BDD-RA-Role-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-RA-Role-{time} members\" dialog"))));
      await session.step(64, "Then the \"BDD-RA-Role-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time} members")));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
      await session.step(66, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Removing the assignment takes the role away (Roles-12)", async () => {
      await session.step(69, "When user clicks on MANAGE button in \"Assigned to\" section in context panel", () => clickOn(page, el("MANAGE button in \"Assigned to\" section in context panel")));
      await session.step(70, "Then \"opavlenko{time}r\" membership row should be visible", () => shouldBe(page, el(session.text("\"opavlenko{time}r\" membership row")), "visible"));
      await session.step(71, "When user clicks on remove button of \"opavlenko{time}r\" membership row", () => clickOn(page, el(session.text("remove button of \"opavlenko{time}r\" membership row"))));
      await session.step(72, "Then \"opavlenko{time}r\" membership row should be absent", () => shouldBe(page, el(session.text("\"opavlenko{time}r\" membership row")), "absent"));
      await session.step(73, "When user clicks on SAVE button in \"BDD-RA-Role-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-RA-Role-{time} members\" dialog"))));
      await session.step(74, "Then the \"BDD-RA-Role-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time} members")));
      await session.step(75, "And \"opavlenko{time}r\" should not be a member of \"BDD-RA-Role-{time}\" on the server", () => notMemberOnServer(page, session.text("opavlenko{time}r"), session.text("BDD-RA-Role-{time}")));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
      await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A new role has no global permissions in its pane (Roles-14)", async () => {
      await session.step(83, "When user expands \"Global Permissions\" section in context panel", () => expand(page, el("\"Global Permissions\" section in context panel")));
      await session.step(84, "Then \"Global Permissions\" section in context panel should contain text \"No global permissions\"", () => shouldContainText(page, el("\"Global Permissions\" section in context panel"), "No global permissions"));
    }, {knownFailure: true});
    await run.scenario("Global Permissions grants the role a permission (Roles-14)", async () => {
      await session.step(87, "When user expands \"Global Permissions\" section in context panel", () => expand(page, el("\"Global Permissions\" section in context panel")));
      await session.step(88, "And user clicks on MANAGE button in \"Global Permissions\" section in context panel", () => clickOn(page, el("MANAGE button in \"Global Permissions\" section in context panel")));
      await session.step(89, "Then \"BDD-RA-Role-{time}: Global Permissions\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-RA-Role-{time}: Global Permissions\" dialog")), "visible"));
      await session.step(90, "When user expands \"Browse\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => expand(page, el(session.text("\"Browse\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(91, "Then \"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog should be unchecked", () => shouldBe(page, el(session.text("\"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog")), "unchecked"));
      await session.step(92, "When user checks \"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => check(page, el(session.text("\"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(93, "And user clicks on SAVE button in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(94, "Then the \"BDD-RA-Role-{time}: Global Permissions\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time}: Global Permissions")));
      await session.step(95, "When user clicks on MANAGE button in \"Global Permissions\" section in context panel", () => clickOn(page, el("MANAGE button in \"Global Permissions\" section in context panel")));
      await session.step(96, "And user expands \"Browse\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => expand(page, el(session.text("\"Browse\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(97, "Then \"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog should be checked", () => shouldBe(page, el(session.text("\"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog")), "checked"));
      await session.step(98, "When user clicks on CANCEL button in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(99, "Then the \"BDD-RA-Role-{time}: Global Permissions\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time}: Global Permissions")));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
      await session.step(101, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A role that grants a permission can still be deleted (Roles-14, Roles-15)", async () => {
      await session.step(110, "When user types \"BDD-RA-Role-{time}\" into gallery search", () => typeInto(page, session.text("BDD-RA-Role-{time}"), el("gallery search")));
      await session.step(111, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(112, "When user picks \"Delete\" from the context menu of \"BDD-RA-Role-{time}\" link in gallery", () => pickFromContextMenu(page, "Delete", el(session.text("\"BDD-RA-Role-{time}\" link in gallery"))));
      await session.step(113, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(114, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(115, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(116, "And 0 roles named \"BDD-RA-Role-{time}\" should be on the server", () => rolesOnServer(page, 0, session.text("BDD-RA-Role-{time}")));
      await session.step(117, "And no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true});
    run.finish();
  });
});
