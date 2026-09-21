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
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, galleryCountLower, groupOnServer, groupsOnServer, memberOnServer, noGroupOnServer, notMemberOnServer, rememberGalleryCount, userOnServer, userStatusOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuDoesNotList, menuLists, noBalloons, noErrors, openContextMenu, pickFromContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Managing a user", () => {
  const session = feature(test, "features/users-groups-roles/users-manage.feature", import.meta.url);
  test("Managing a user", {tag: ["@journey", "@serial", "@users", "@realizes:views.users"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And a user \"bddmanaged\" is on the server", () => userOnServer(page, "bddmanaged"));
    await session.step(24, "And a group named \"BDD-UM-Group-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-UM-Group-{time}")));
    await session.step(25, "And no role named \"BDD-UM-Role-{time}\" is on the server", () => noGroupOnServer(page, session.text("BDD-UM-Role-{time}")));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await run.scenario("A role to assign", async () => {
      await session.step(31, "When user clicks on \"Platform > Roles\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Roles\" tree node inside browse tree")));
      await session.step(32, "Then the \"Roles\" view should be current", () => viewIsCurrent(page, "Roles"));
      await session.step(33, "When user clicks on \"New Role...\" button", () => clickOn(page, el("\"New Role...\" button")));
      await session.step(34, "And user types \"BDD-UM-Role-{time}\" into Name input in \"Create New Role\" dialog", () => typeInto(page, session.text("BDD-UM-Role-{time}"), el("Name input in \"Create New Role\" dialog")));
      await session.step(35, "And user clicks on OK button in \"Create New Role\" dialog", () => clickOn(page, el("OK button in \"Create New Role\" dialog")));
      await session.step(36, "Then the \"Create New Role\" dialog should close", () => dialogCloses(page, "Create New Role"));
      await session.step(37, "And 1 role named \"BDD-UM-Role-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-UM-Role-{time}")));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
      await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Groups... adds the user to a group (Users-18)", async () => {
      await session.step(42, "When user clicks on \"Platform > Users\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Users\" tree node inside browse tree")));
      await session.step(43, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(44, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(45, "And user types \"bddmanaged\" into gallery search", () => typeInto(page, "bddmanaged", el("gallery search")));
      await session.step(46, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(47, "When user picks \"Groups...\" from the context menu of \"bddmanaged\" link in gallery", () => pickFromContextMenu(page, "Groups...", el("\"bddmanaged\" link in gallery")));
      await session.step(48, "Then \"bddmanaged groups\" dialog should be visible", () => shouldBe(page, el("\"bddmanaged groups\" dialog"), "visible"));
      await session.step(49, "When user types \"BDD-UM-Group-{time}\" into membership search", () => typeInto(page, session.text("BDD-UM-Group-{time}"), el("membership search")));
      await session.step(50, "And user clicks on add button of \"BDD-UM-Group-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-UM-Group-{time}\" membership candidate"))));
      await session.step(51, "Then \"BDD-UM-Group-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-UM-Group-{time}\" membership row")), "visible"));
      await session.step(52, "When user clicks on SAVE button in \"bddmanaged groups\" dialog", () => clickOn(page, el("SAVE button in \"bddmanaged groups\" dialog")));
      await session.step(53, "Then the \"bddmanaged groups\" dialog should close", () => dialogCloses(page, "bddmanaged groups"));
      await session.step(54, "And \"bddmanaged\" should be a member of \"BDD-UM-Group-{time}\" on the server", () => memberOnServer(page, "bddmanaged", session.text("BDD-UM-Group-{time}")));
      await session.step(55, "When user clicks on \"bddmanaged\" link in gallery", () => clickOn(page, el("\"bddmanaged\" link in gallery")));
      await session.step(56, "Then the context panel should show \"bddmanaged\"", () => contextPanelShows(page, "bddmanaged"));
      await session.step(57, "When user expands \"Member of\" section in context panel", () => expand(page, el("\"Member of\" section in context panel")));
      await session.step(58, "Then \"Member of\" section in context panel should contain text \"BDD-UM-Group-{time}\"", () => shouldContainText(page, el("\"Member of\" section in context panel"), session.text("BDD-UM-Group-{time}")));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
      await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Groups... takes the user out of the group again (Users-18)", async () => {
      await session.step(63, "When user picks \"Groups...\" from the context menu of \"bddmanaged\" link in gallery", () => pickFromContextMenu(page, "Groups...", el("\"bddmanaged\" link in gallery")));
      await session.step(64, "Then \"BDD-UM-Group-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-UM-Group-{time}\" membership row")), "visible"));
      await session.step(65, "When user clicks on remove button of \"BDD-UM-Group-{time}\" membership row", () => clickOn(page, el(session.text("remove button of \"BDD-UM-Group-{time}\" membership row"))));
      await session.step(66, "Then \"BDD-UM-Group-{time}\" membership row should be absent", () => shouldBe(page, el(session.text("\"BDD-UM-Group-{time}\" membership row")), "absent"));
      await session.step(67, "When user clicks on SAVE button in \"bddmanaged groups\" dialog", () => clickOn(page, el("SAVE button in \"bddmanaged groups\" dialog")));
      await session.step(68, "Then the \"bddmanaged groups\" dialog should close", () => dialogCloses(page, "bddmanaged groups"));
      await session.step(69, "And \"bddmanaged\" should not be a member of \"BDD-UM-Group-{time}\" on the server", () => notMemberOnServer(page, "bddmanaged", session.text("BDD-UM-Group-{time}")));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
      await session.step(71, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Roles... gives the user a role (Users-19)", async () => {
      await session.step(74, "When user picks \"Roles...\" from the context menu of \"bddmanaged\" link in gallery", () => pickFromContextMenu(page, "Roles...", el("\"bddmanaged\" link in gallery")));
      await session.step(75, "Then \"bddmanaged roles\" dialog should be visible", () => shouldBe(page, el("\"bddmanaged roles\" dialog"), "visible"));
      await session.step(76, "When user types \"BDD-UM-Role-{time}\" into membership search", () => typeInto(page, session.text("BDD-UM-Role-{time}"), el("membership search")));
      await session.step(77, "And user clicks on add button of \"BDD-UM-Role-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-UM-Role-{time}\" membership candidate"))));
      await session.step(78, "Then \"BDD-UM-Role-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-UM-Role-{time}\" membership row")), "visible"));
      await session.step(79, "When user clicks on SAVE button in \"bddmanaged roles\" dialog", () => clickOn(page, el("SAVE button in \"bddmanaged roles\" dialog")));
      await session.step(80, "Then the \"bddmanaged roles\" dialog should close", () => dialogCloses(page, "bddmanaged roles"));
      await session.step(81, "And \"bddmanaged\" should be a member of \"BDD-UM-Role-{time}\" on the server", () => memberOnServer(page, "bddmanaged", session.text("BDD-UM-Role-{time}")));
      await session.step(82, "When user clicks on \"bddmanaged\" link in gallery", () => clickOn(page, el("\"bddmanaged\" link in gallery")));
      await session.step(83, "And user expands \"Roles\" section in context panel", () => expand(page, el("\"Roles\" section in context panel")));
      await session.step(84, "Then \"Roles\" section in context panel should contain text \"BDD-UM-Role-{time}\"", () => shouldContainText(page, el("\"Roles\" section in context panel"), session.text("BDD-UM-Role-{time}")));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
      await session.step(86, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Roles... takes the role away again (Users-19)", async () => {
      await session.step(89, "When user picks \"Roles...\" from the context menu of \"bddmanaged\" link in gallery", () => pickFromContextMenu(page, "Roles...", el("\"bddmanaged\" link in gallery")));
      await session.step(90, "Then \"BDD-UM-Role-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-UM-Role-{time}\" membership row")), "visible"));
      await session.step(91, "When user clicks on remove button of \"BDD-UM-Role-{time}\" membership row", () => clickOn(page, el(session.text("remove button of \"BDD-UM-Role-{time}\" membership row"))));
      await session.step(92, "Then \"BDD-UM-Role-{time}\" membership row should be absent", () => shouldBe(page, el(session.text("\"BDD-UM-Role-{time}\" membership row")), "absent"));
      await session.step(93, "When user clicks on SAVE button in \"bddmanaged roles\" dialog", () => clickOn(page, el("SAVE button in \"bddmanaged roles\" dialog")));
      await session.step(94, "Then the \"bddmanaged roles\" dialog should close", () => dialogCloses(page, "bddmanaged roles"));
      await session.step(95, "And \"bddmanaged\" should not be a member of \"BDD-UM-Role-{time}\" on the server", () => notMemberOnServer(page, "bddmanaged", session.text("BDD-UM-Role-{time}")));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
      await session.step(97, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Disable... disables the user (Users-20)", async () => {
      await session.step(100, "Then the user \"bddmanaged\" should be active on the server", () => userStatusOnServer(page, "bddmanaged", "active"));
      await session.step(101, "When user picks \"Disable...\" from the context menu of \"bddmanaged\" link in gallery", () => pickFromContextMenu(page, "Disable...", el("\"bddmanaged\" link in gallery")));
      await session.step(102, "Then \"Disable user\" dialog should be visible", () => shouldBe(page, el("\"Disable user\" dialog"), "visible"));
      await session.step(103, "When user clicks on DISABLE button in \"Disable user\" dialog", () => clickOn(page, el("DISABLE button in \"Disable user\" dialog")));
      await session.step(104, "Then the \"Disable user\" dialog should close", () => dialogCloses(page, "Disable user"));
      await session.step(105, "And the user \"bddmanaged\" should be disabled on the server", () => userStatusOnServer(page, "bddmanaged", "disabled"));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
      await session.step(107, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A disabled user can be enabled again (Users-20)", async () => {
      await session.step(110, "When user opens the context menu of \"bddmanaged\" link in gallery", () => openContextMenu(page, el("\"bddmanaged\" link in gallery")));
      await session.step(111, "Then the open menu should list \"Enable\"", () => menuLists(page, "Enable"));
      await session.step(112, "And the open menu should not list \"Disable...\"", () => menuDoesNotList(page, "Disable..."));
      await session.step(113, "When user picks \"Enable\" from the open menu", () => pickFromOpenMenu(page, "Enable"));
      await session.step(114, "Then \"Enable user\" dialog should be visible", () => shouldBe(page, el("\"Enable user\" dialog"), "visible"));
      await session.step(115, "When user clicks on ENABLE button in \"Enable user\" dialog", () => clickOn(page, el("ENABLE button in \"Enable user\" dialog")));
      await session.step(116, "Then the \"Enable user\" dialog should close", () => dialogCloses(page, "Enable user"));
      await session.step(117, "And the user \"bddmanaged\" should be active on the server", () => userStatusOnServer(page, "bddmanaged", "active"));
      await session.step(118, "And no errors should have been logged", () => noErrors(page));
      await session.step(119, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A user is added to favorites and removed again (Users-21)", async () => {
      await session.step(122, "Then \"My stuff > Favorites > bddmanaged\" tree node inside browse tree should be absent", () => shouldBe(page, el("\"My stuff > Favorites > bddmanaged\" tree node inside browse tree"), "absent"));
      await session.step(123, "When user picks \"Add to favorites\" from the context menu of \"bddmanaged\" link in gallery", () => pickFromContextMenu(page, "Add to favorites", el("\"bddmanaged\" link in gallery")));
      await session.step(124, "Then \"My stuff > Favorites > bddmanaged\" tree node inside browse tree should be present", () => shouldBe(page, el("\"My stuff > Favorites > bddmanaged\" tree node inside browse tree"), "present"));
      await session.step(125, "When user opens the context menu of \"bddmanaged\" link in gallery", () => openContextMenu(page, el("\"bddmanaged\" link in gallery")));
      await session.step(126, "Then the open menu should list \"Remove from favorites\"", () => menuLists(page, "Remove from favorites"));
      await session.step(127, "When user picks \"Remove from favorites\" from the open menu", () => pickFromOpenMenu(page, "Remove from favorites"));
      await session.step(128, "Then \"My stuff > Favorites > bddmanaged\" tree node inside browse tree should be absent", () => shouldBe(page, el("\"My stuff > Favorites > bddmanaged\" tree node inside browse tree"), "absent"));
      await session.step(129, "When user opens the context menu of \"bddmanaged\" link in gallery", () => openContextMenu(page, el("\"bddmanaged\" link in gallery")));
      await session.step(130, "Then the open menu should list \"Add to favorites\"", () => menuLists(page, "Add to favorites"));
      await session.step(131, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(132, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(133, "Then no errors should have been logged", () => noErrors(page));
      await session.step(134, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
