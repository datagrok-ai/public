/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/users-manage.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.users]
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
import {clearField, clickOn, expand, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, galleryCountLower, groupOnServer, memberOnServer, newUserOnServer, noRoleOnServer, notMemberOnServer, rememberGalleryCount, rolesOnServer, userStatusOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuDoesNotList, menuLists, noBalloons, noErrors, openContextMenu, pickFromContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Managing a user", () => {
  const session = feature(test, "features/users-groups-roles/users-manage.feature", import.meta.url);
  test("Managing a user", {tag: ["@journey", "@serial", "@users", "@realizes:views.users"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And a new user \"opavlenko{time}m\" with email \"opavlenko+{time}m@datagrok.ai\" is on the server", () => newUserOnServer(page, session.text("opavlenko{time}m"), session.text("opavlenko+{time}m@datagrok.ai")));
    await session.step(23, "And a group named \"BDD-UM-Group-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-UM-Group-{time}")));
    await session.step(24, "And no role named \"BDD-UM-Role-{time}\" is on the server", () => noRoleOnServer(page, session.text("BDD-UM-Role-{time}")));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(27, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await run.scenario("A role to assign", async () => {
      await session.step(30, "When user clicks on \"Platform > Roles\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Roles\" tree node inside browse tree")));
      await session.step(31, "Then the \"Roles\" view should be current", () => viewIsCurrent(page, "Roles"));
      await session.step(32, "When user clicks on \"New Role...\" button", () => clickOn(page, el("\"New Role...\" button")));
      await session.step(33, "And user types \"BDD-UM-Role-{time}\" into Name input in \"Create New Role\" dialog", () => typeInto(page, session.text("BDD-UM-Role-{time}"), el("Name input in \"Create New Role\" dialog")));
      await session.step(34, "And user clicks on OK button in \"Create New Role\" dialog", () => clickOn(page, el("OK button in \"Create New Role\" dialog")));
      await session.step(35, "Then the \"Create New Role\" dialog should close", () => dialogCloses(page, "Create New Role"));
      await session.step(36, "And 1 role named \"BDD-UM-Role-{time}\" should be on the server", () => rolesOnServer(page, 1, session.text("BDD-UM-Role-{time}")));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
      await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Groups... adds the user to a group (Users-18)", async () => {
      await session.step(41, "When user clicks on \"Platform > Users\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Users\" tree node inside browse tree")));
      await session.step(42, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(43, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(44, "And user types \"opavlenko{time}m\" into gallery search", () => typeInto(page, session.text("opavlenko{time}m"), el("gallery search")));
      await session.step(45, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(46, "When user picks \"Groups...\" from the context menu of \"opavlenko{time}m\" link in gallery", () => pickFromContextMenu(page, "Groups...", el(session.text("\"opavlenko{time}m\" link in gallery"))));
      await session.step(47, "Then \"opavlenko{time}m groups\" dialog should be visible", () => shouldBe(page, el(session.text("\"opavlenko{time}m groups\" dialog")), "visible"));
      await session.step(48, "When user types \"BDD-UM-Group-{time}\" into membership search", () => typeInto(page, session.text("BDD-UM-Group-{time}"), el("membership search")));
      await session.step(49, "And user clicks on add button of \"BDD-UM-Group-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-UM-Group-{time}\" membership candidate"))));
      await session.step(50, "Then \"BDD-UM-Group-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-UM-Group-{time}\" membership row")), "visible"));
      await session.step(51, "When user clicks on SAVE button in \"opavlenko{time}m groups\" dialog", () => clickOn(page, el(session.text("SAVE button in \"opavlenko{time}m groups\" dialog"))));
      await session.step(52, "Then the \"opavlenko{time}m groups\" dialog should close", () => dialogCloses(page, session.text("opavlenko{time}m groups")));
      await session.step(53, "And \"opavlenko{time}m\" should be a member of \"BDD-UM-Group-{time}\" on the server", () => memberOnServer(page, session.text("opavlenko{time}m"), session.text("BDD-UM-Group-{time}")));
      await session.step(54, "When user clicks on \"opavlenko{time}m\" link in gallery", () => clickOn(page, el(session.text("\"opavlenko{time}m\" link in gallery"))));
      await session.step(55, "Then the context panel should show \"opavlenko{time}m\"", () => contextPanelShows(page, session.text("opavlenko{time}m")));
      await session.step(56, "When user expands \"Member of\" section in context panel", () => expand(page, el("\"Member of\" section in context panel")));
      await session.step(57, "Then \"Member of\" section in context panel should contain text \"BDD-UM-Group-{time}\"", () => shouldContainText(page, el("\"Member of\" section in context panel"), session.text("BDD-UM-Group-{time}")));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
      await session.step(59, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Groups... takes the user out of the group again (Users-18)", async () => {
      await session.step(62, "When user picks \"Groups...\" from the context menu of \"opavlenko{time}m\" link in gallery", () => pickFromContextMenu(page, "Groups...", el(session.text("\"opavlenko{time}m\" link in gallery"))));
      await session.step(63, "Then \"BDD-UM-Group-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-UM-Group-{time}\" membership row")), "visible"));
      await session.step(64, "When user clicks on remove button of \"BDD-UM-Group-{time}\" membership row", () => clickOn(page, el(session.text("remove button of \"BDD-UM-Group-{time}\" membership row"))));
      await session.step(65, "Then \"BDD-UM-Group-{time}\" membership row should be absent", () => shouldBe(page, el(session.text("\"BDD-UM-Group-{time}\" membership row")), "absent"));
      await session.step(66, "When user clicks on SAVE button in \"opavlenko{time}m groups\" dialog", () => clickOn(page, el(session.text("SAVE button in \"opavlenko{time}m groups\" dialog"))));
      await session.step(67, "Then the \"opavlenko{time}m groups\" dialog should close", () => dialogCloses(page, session.text("opavlenko{time}m groups")));
      await session.step(68, "And \"opavlenko{time}m\" should not be a member of \"BDD-UM-Group-{time}\" on the server", () => notMemberOnServer(page, session.text("opavlenko{time}m"), session.text("BDD-UM-Group-{time}")));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
      await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Roles... gives the user a role (Users-19)", async () => {
      await session.step(73, "When user picks \"Roles...\" from the context menu of \"opavlenko{time}m\" link in gallery", () => pickFromContextMenu(page, "Roles...", el(session.text("\"opavlenko{time}m\" link in gallery"))));
      await session.step(74, "Then \"opavlenko{time}m roles\" dialog should be visible", () => shouldBe(page, el(session.text("\"opavlenko{time}m roles\" dialog")), "visible"));
      await session.step(75, "When user types \"BDD-UM-Role-{time}\" into membership search", () => typeInto(page, session.text("BDD-UM-Role-{time}"), el("membership search")));
      await session.step(76, "And user clicks on add button of \"BDD-UM-Role-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-UM-Role-{time}\" membership candidate"))));
      await session.step(77, "Then \"BDD-UM-Role-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-UM-Role-{time}\" membership row")), "visible"));
      await session.step(78, "When user clicks on SAVE button in \"opavlenko{time}m roles\" dialog", () => clickOn(page, el(session.text("SAVE button in \"opavlenko{time}m roles\" dialog"))));
      await session.step(79, "Then the \"opavlenko{time}m roles\" dialog should close", () => dialogCloses(page, session.text("opavlenko{time}m roles")));
      await session.step(80, "And \"opavlenko{time}m\" should be a member of \"BDD-UM-Role-{time}\" on the server", () => memberOnServer(page, session.text("opavlenko{time}m"), session.text("BDD-UM-Role-{time}")));
      await session.step(81, "When user clicks on \"opavlenko{time}m\" link in gallery", () => clickOn(page, el(session.text("\"opavlenko{time}m\" link in gallery"))));
      await session.step(82, "And user expands \"Roles\" section in context panel", () => expand(page, el("\"Roles\" section in context panel")));
      await session.step(83, "Then \"Roles\" section in context panel should contain text \"BDD-UM-Role-{time}\"", () => shouldContainText(page, el("\"Roles\" section in context panel"), session.text("BDD-UM-Role-{time}")));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
      await session.step(85, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Roles... takes the role away again (Users-19)", async () => {
      await session.step(88, "When user picks \"Roles...\" from the context menu of \"opavlenko{time}m\" link in gallery", () => pickFromContextMenu(page, "Roles...", el(session.text("\"opavlenko{time}m\" link in gallery"))));
      await session.step(89, "Then \"BDD-UM-Role-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-UM-Role-{time}\" membership row")), "visible"));
      await session.step(90, "When user clicks on remove button of \"BDD-UM-Role-{time}\" membership row", () => clickOn(page, el(session.text("remove button of \"BDD-UM-Role-{time}\" membership row"))));
      await session.step(91, "Then \"BDD-UM-Role-{time}\" membership row should be absent", () => shouldBe(page, el(session.text("\"BDD-UM-Role-{time}\" membership row")), "absent"));
      await session.step(92, "When user clicks on SAVE button in \"opavlenko{time}m roles\" dialog", () => clickOn(page, el(session.text("SAVE button in \"opavlenko{time}m roles\" dialog"))));
      await session.step(93, "Then the \"opavlenko{time}m roles\" dialog should close", () => dialogCloses(page, session.text("opavlenko{time}m roles")));
      await session.step(94, "And \"opavlenko{time}m\" should not be a member of \"BDD-UM-Role-{time}\" on the server", () => notMemberOnServer(page, session.text("opavlenko{time}m"), session.text("BDD-UM-Role-{time}")));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
      await session.step(96, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Disable... disables the user (Users-20)", async () => {
      await session.step(99, "Then the user \"opavlenko{time}m\" should be active on the server", () => userStatusOnServer(page, session.text("opavlenko{time}m"), "active"));
      await session.step(100, "When user picks \"Disable...\" from the context menu of \"opavlenko{time}m\" link in gallery", () => pickFromContextMenu(page, "Disable...", el(session.text("\"opavlenko{time}m\" link in gallery"))));
      await session.step(101, "Then \"Disable user\" dialog should be visible", () => shouldBe(page, el("\"Disable user\" dialog"), "visible"));
      await session.step(102, "When user clicks on DISABLE button in \"Disable user\" dialog", () => clickOn(page, el("DISABLE button in \"Disable user\" dialog")));
      await session.step(103, "Then the \"Disable user\" dialog should close", () => dialogCloses(page, "Disable user"));
      await session.step(104, "And the user \"opavlenko{time}m\" should be disabled on the server", () => userStatusOnServer(page, session.text("opavlenko{time}m"), "disabled"));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
      await session.step(106, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A disabled user can be enabled again (Users-20)", async () => {
      await session.step(109, "When user opens the context menu of \"opavlenko{time}m\" link in gallery", () => openContextMenu(page, el(session.text("\"opavlenko{time}m\" link in gallery"))));
      await session.step(110, "Then the open menu should list \"Enable\"", () => menuLists(page, "Enable"));
      await session.step(111, "And the open menu should not list \"Disable...\"", () => menuDoesNotList(page, "Disable..."));
      await session.step(112, "When user picks \"Enable\" from the open menu", () => pickFromOpenMenu(page, "Enable"));
      await session.step(113, "Then \"Enable user\" dialog should be visible", () => shouldBe(page, el("\"Enable user\" dialog"), "visible"));
      await session.step(114, "When user clicks on ENABLE button in \"Enable user\" dialog", () => clickOn(page, el("ENABLE button in \"Enable user\" dialog")));
      await session.step(115, "Then the \"Enable user\" dialog should close", () => dialogCloses(page, "Enable user"));
      await session.step(116, "And the user \"opavlenko{time}m\" should be active on the server", () => userStatusOnServer(page, session.text("opavlenko{time}m"), "active"));
      await session.step(117, "And no errors should have been logged", () => noErrors(page));
      await session.step(118, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A user is added to favorites and removed again (Users-21)", async () => {
      await session.step(121, "Then \"My stuff > Favorites > opavlenko{time}m\" tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("\"My stuff > Favorites > opavlenko{time}m\" tree node inside browse tree")), "absent"));
      await session.step(122, "When user picks \"Add to favorites\" from the context menu of \"opavlenko{time}m\" link in gallery", () => pickFromContextMenu(page, "Add to favorites", el(session.text("\"opavlenko{time}m\" link in gallery"))));
      await session.step(123, "Then \"My stuff > Favorites > opavlenko{time}m\" tree node inside browse tree should be present", () => shouldBe(page, el(session.text("\"My stuff > Favorites > opavlenko{time}m\" tree node inside browse tree")), "present"));
      await session.step(124, "When user opens the context menu of \"opavlenko{time}m\" link in gallery", () => openContextMenu(page, el(session.text("\"opavlenko{time}m\" link in gallery"))));
      await session.step(125, "Then the open menu should list \"Remove from favorites\"", () => menuLists(page, "Remove from favorites"));
      await session.step(126, "When user picks \"Remove from favorites\" from the open menu", () => pickFromOpenMenu(page, "Remove from favorites"));
      await session.step(127, "Then \"My stuff > Favorites > opavlenko{time}m\" tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("\"My stuff > Favorites > opavlenko{time}m\" tree node inside browse tree")), "absent"));
      await session.step(128, "When user opens the context menu of \"opavlenko{time}m\" link in gallery", () => openContextMenu(page, el(session.text("\"opavlenko{time}m\" link in gallery"))));
      await session.step(129, "Then the open menu should list \"Add to favorites\"", () => menuLists(page, "Add to favorites"));
      await session.step(130, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(131, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(132, "Then no errors should have been logged", () => noErrors(page));
      await session.step(133, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
