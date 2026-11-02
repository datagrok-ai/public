/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/roles-assignment.feature
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
import {check, clickOn, expand, shouldBe, shouldContainText, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {adminMemberOnServer, browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, galleryCountHigher, galleryCountLower, groupsOnServer, noGroupOnServer, notMemberOnServer, plainMemberOnServer, rememberGalleryCount, userOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Who holds a role, and what it grants", () => {
  const session = feature(test, "features/users-groups-roles/roles-assignment.feature", import.meta.url);
  test("Who holds a role, and what it grants", {tag: ["@journey", "@serial", "@roles", "@realizes:views.roles"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And a user \"bddmanaged\" is on the server", () => userOnServer(page, "bddmanaged"));
    await session.step(27, "And no role named \"BDD-RA-Role-{time}\" is on the server", () => noGroupOnServer(page, session.text("BDD-RA-Role-{time}")));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(30, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(31, "And user clicks on \"Platform > Roles\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Roles\" tree node inside browse tree")));
    await run.scenario("A role to assign", async () => {
      await session.step(34, "Then the \"Roles\" view should be current", () => viewIsCurrent(page, "Roles"));
      await session.step(35, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(36, "And user clicks on \"New Role...\" button", () => clickOn(page, el("\"New Role...\" button")));
      await session.step(37, "And user types \"BDD-RA-Role-{time}\" into Name input in \"Create New Role\" dialog", () => typeInto(page, session.text("BDD-RA-Role-{time}"), el("Name input in \"Create New Role\" dialog")));
      await session.step(38, "And user clicks on OK button in \"Create New Role\" dialog", () => clickOn(page, el("OK button in \"Create New Role\" dialog")));
      await session.step(39, "Then the \"Create New Role\" dialog should close", () => dialogCloses(page, "Create New Role"));
      await session.step(40, "And 1 role named \"BDD-RA-Role-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-RA-Role-{time}")));
      await session.step(41, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(42, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
      await session.step(43, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(44, "And user types \"BDD-RA-Role-{time}\" into gallery search", () => typeInto(page, session.text("BDD-RA-Role-{time}"), el("gallery search")));
      await session.step(45, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(46, "When user clicks on \"BDD-RA-Role-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-RA-Role-{time}\" link in gallery"))));
      await session.step(47, "Then the context panel should show \"BDD-RA-Role-{time}\"", () => contextPanelShows(page, session.text("BDD-RA-Role-{time}")));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
      await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("MANAGE assigns the role to a user (Roles-11)", async () => {
      await session.step(52, "When user clicks on MANAGE button in \"Assigned to\" section in context panel", () => clickOn(page, el("MANAGE button in \"Assigned to\" section in context panel")));
      await session.step(53, "Then \"BDD-RA-Role-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-RA-Role-{time} members\" dialog")), "visible"));
      await session.step(54, "When user types \"bddmanaged\" into membership search", () => typeInto(page, "bddmanaged", el("membership search")));
      await session.step(55, "And user clicks on add button of \"bddmanaged\" membership candidate", () => clickOn(page, el("add button of \"bddmanaged\" membership candidate")));
      await session.step(56, "Then \"bddmanaged\" membership row should be visible", () => shouldBe(page, el("\"bddmanaged\" membership row"), "visible"));
      await session.step(57, "And checkbox label of \"bddmanaged\" membership row should have text \"Can assign\"", () => shouldHaveText(page, el("checkbox label of \"bddmanaged\" membership row"), "Can assign"));
      await session.step(58, "When user clicks on SAVE button in \"BDD-RA-Role-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-RA-Role-{time} members\" dialog"))));
      await session.step(59, "Then the \"BDD-RA-Role-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time} members")));
      await session.step(60, "And \"bddmanaged\" should be a plain member of \"BDD-RA-Role-{time}\" on the server", () => plainMemberOnServer(page, "bddmanaged", session.text("BDD-RA-Role-{time}")));
      await session.step(61, "And \"Assigned to\" section in context panel should contain text \"bddmanaged\"", () => shouldContainText(page, el("\"Assigned to\" section in context panel"), "bddmanaged"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
      await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Can assign is saved for the assignee (Roles-13)", async () => {
      await session.step(66, "When user clicks on MANAGE button in \"Assigned to\" section in context panel", () => clickOn(page, el("MANAGE button in \"Assigned to\" section in context panel")));
      await session.step(67, "Then checkbox of \"bddmanaged\" membership row should be unchecked", () => shouldBe(page, el("checkbox of \"bddmanaged\" membership row"), "unchecked"));
      await session.step(68, "When user checks checkbox of \"bddmanaged\" membership row", () => check(page, el("checkbox of \"bddmanaged\" membership row")));
      await session.step(69, "And user clicks on SAVE button in \"BDD-RA-Role-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-RA-Role-{time} members\" dialog"))));
      await session.step(70, "Then the \"BDD-RA-Role-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time} members")));
      await session.step(71, "And \"bddmanaged\" should be an admin member of \"BDD-RA-Role-{time}\" on the server", () => adminMemberOnServer(page, "bddmanaged", session.text("BDD-RA-Role-{time}")));
      await session.step(72, "When user clicks on MANAGE button in \"Assigned to\" section in context panel", () => clickOn(page, el("MANAGE button in \"Assigned to\" section in context panel")));
      await session.step(73, "Then checkbox of \"bddmanaged\" membership row should be checked", () => shouldBe(page, el("checkbox of \"bddmanaged\" membership row"), "checked"));
      await session.step(74, "When user clicks on CANCEL button in \"BDD-RA-Role-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-RA-Role-{time} members\" dialog"))));
      await session.step(75, "Then the \"BDD-RA-Role-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time} members")));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
      await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Removing the assignment takes the role away (Roles-12)", async () => {
      await session.step(80, "When user clicks on MANAGE button in \"Assigned to\" section in context panel", () => clickOn(page, el("MANAGE button in \"Assigned to\" section in context panel")));
      await session.step(81, "Then \"bddmanaged\" membership row should be visible", () => shouldBe(page, el("\"bddmanaged\" membership row"), "visible"));
      await session.step(82, "When user clicks on remove button of \"bddmanaged\" membership row", () => clickOn(page, el("remove button of \"bddmanaged\" membership row")));
      await session.step(83, "Then \"bddmanaged\" membership row should be absent", () => shouldBe(page, el("\"bddmanaged\" membership row"), "absent"));
      await session.step(84, "When user clicks on SAVE button in \"BDD-RA-Role-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-RA-Role-{time} members\" dialog"))));
      await session.step(85, "Then the \"BDD-RA-Role-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time} members")));
      await session.step(86, "And \"bddmanaged\" should not be a member of \"BDD-RA-Role-{time}\" on the server", () => notMemberOnServer(page, "bddmanaged", session.text("BDD-RA-Role-{time}")));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
      await session.step(88, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A new role has no global permissions in its pane (Roles-14)", async () => {
      await session.step(93, "When user expands \"Global Permissions\" section in context panel", () => expand(page, el("\"Global Permissions\" section in context panel")));
      await session.step(94, "Then \"Global Permissions\" section in context panel should contain text \"No global permissions\"", () => shouldContainText(page, el("\"Global Permissions\" section in context panel"), "No global permissions"));
      await session.step(95, "When user clicks on MANAGE button in \"Global Permissions\" section in context panel", () => clickOn(page, el("MANAGE button in \"Global Permissions\" section in context panel")));
      await session.step(96, "Then \"BDD-RA-Role-{time}: Global Permissions\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-RA-Role-{time}: Global Permissions\" dialog")), "visible"));
      await session.step(97, "When user expands \"Browse\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => expand(page, el(session.text("\"Browse\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(98, "And user checks \"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => check(page, el(session.text("\"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(99, "Then \"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog should be checked", () => shouldBe(page, el(session.text("\"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog")), "checked"));
      await session.step(100, "When user clicks on CANCEL button in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(101, "Then the \"BDD-RA-Role-{time}: Global Permissions\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time}: Global Permissions")));
      await session.step(102, "When user clicks on MANAGE button in \"Global Permissions\" section in context panel", () => clickOn(page, el("MANAGE button in \"Global Permissions\" section in context panel")));
      await session.step(103, "And user expands \"Browse\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => expand(page, el(session.text("\"Browse\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(104, "Then \"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog should be unchecked", () => shouldBe(page, el(session.text("\"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog")), "unchecked"));
      await session.step(105, "When user clicks on CANCEL button in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(106, "Then the \"BDD-RA-Role-{time}: Global Permissions\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time}: Global Permissions")));
      await session.step(107, "And \"Global Permissions\" section in context panel should contain text \"No global permissions\"", () => shouldContainText(page, el("\"Global Permissions\" section in context panel"), "No global permissions"));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
      await session.step(109, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Global Permissions grants the role a permission (Roles-14)", async () => {
      await session.step(112, "When user expands \"Global Permissions\" section in context panel", () => expand(page, el("\"Global Permissions\" section in context panel")));
      await session.step(113, "And user clicks on MANAGE button in \"Global Permissions\" section in context panel", () => clickOn(page, el("MANAGE button in \"Global Permissions\" section in context panel")));
      await session.step(114, "Then \"BDD-RA-Role-{time}: Global Permissions\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-RA-Role-{time}: Global Permissions\" dialog")), "visible"));
      await session.step(115, "When user expands \"Browse\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => expand(page, el(session.text("\"Browse\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(116, "Then \"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog should be unchecked", () => shouldBe(page, el(session.text("\"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog")), "unchecked"));
      await session.step(117, "When user checks \"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => check(page, el(session.text("\"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(118, "And user clicks on SAVE button in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(119, "Then the \"BDD-RA-Role-{time}: Global Permissions\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time}: Global Permissions")));
      await session.step(120, "When user clicks on MANAGE button in \"Global Permissions\" section in context panel", () => clickOn(page, el("MANAGE button in \"Global Permissions\" section in context panel")));
      await session.step(121, "And user expands \"Browse\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => expand(page, el(session.text("\"Browse\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(122, "Then \"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog should be checked", () => shouldBe(page, el(session.text("\"Browse > Browse Apps\" tree node in \"BDD-RA-Role-{time}: Global Permissions\" dialog")), "checked"));
      await session.step(123, "When user clicks on CANCEL button in \"BDD-RA-Role-{time}: Global Permissions\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-RA-Role-{time}: Global Permissions\" dialog"))));
      await session.step(124, "Then the \"BDD-RA-Role-{time}: Global Permissions\" dialog should close", () => dialogCloses(page, session.text("BDD-RA-Role-{time}: Global Permissions")));
      await session.step(125, "And no errors should have been logged", () => noErrors(page));
      await session.step(126, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
