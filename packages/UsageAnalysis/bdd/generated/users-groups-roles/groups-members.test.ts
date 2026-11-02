/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/groups-members.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.groups]
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
import {check, clearField, clickOn, expand, shouldBe, shouldContainText, shouldHaveText, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {adminMemberOnServer, browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, galleryCountLower, groupOnServer, memberOnServer, notMemberOnServer, plainMemberOnServer, rememberGalleryCount, userOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A group's members", () => {
  const session = feature(test, "features/users-groups-roles/groups-members.feature", import.meta.url);
  test("A group's members", {tag: ["@journey", "@serial", "@groups", "@realizes:views.groups", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And a user \"bddmanaged\" is on the server", () => userOnServer(page, "bddmanaged"));
    await session.step(31, "And a group named \"BDD-GM-Group-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Group-{time}")));
    await session.step(32, "And a group named \"BDD-GM-Child-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Child-{time}")));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(34, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(35, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(36, "And user clicks on \"Platform > Groups\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Groups\" tree node inside browse tree")));
    await run.scenario("MANAGE adds a user to the group (Groups-11)", async () => {
      await session.step(39, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(40, "And user types \"BDD-GM-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Group-{time}"), el("gallery search")));
      await session.step(41, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(42, "When user clicks on \"BDD-GM-Group-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Group-{time}\" link in gallery"))));
      await session.step(43, "Then the context panel should show \"BDD-GM-Group-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Group-{time}")));
      await session.step(44, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(45, "Then \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(46, "When user types \"bddmanaged\" into membership search", () => typeInto(page, "bddmanaged", el("membership search")));
      await session.step(47, "And user clicks on add button of \"bddmanaged\" membership candidate", () => clickOn(page, el("add button of \"bddmanaged\" membership candidate")));
      await session.step(48, "Then \"bddmanaged\" membership row should be visible", () => shouldBe(page, el("\"bddmanaged\" membership row"), "visible"));
      await session.step(49, "When user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(50, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(51, "And \"bddmanaged\" should be a plain member of \"BDD-GM-Group-{time}\" on the server", () => plainMemberOnServer(page, "bddmanaged", session.text("BDD-GM-Group-{time}")));
      await session.step(52, "And \"Members\" section in context panel should contain text \"bddmanaged\"", () => shouldContainText(page, el("\"Members\" section in context panel"), "bddmanaged"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Admin box makes the member an admin (Groups-13)", async () => {
      await session.step(57, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(58, "Then checkbox label of \"bddmanaged\" membership row should have text \"Admin\"", () => shouldHaveText(page, el("checkbox label of \"bddmanaged\" membership row"), "Admin"));
      await session.step(59, "And checkbox of \"bddmanaged\" membership row should be unchecked", () => shouldBe(page, el("checkbox of \"bddmanaged\" membership row"), "unchecked"));
      await session.step(60, "When user checks checkbox of \"bddmanaged\" membership row", () => check(page, el("checkbox of \"bddmanaged\" membership row")));
      await session.step(61, "And user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(62, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(63, "And \"bddmanaged\" should be an admin member of \"BDD-GM-Group-{time}\" on the server", () => adminMemberOnServer(page, "bddmanaged", session.text("BDD-GM-Group-{time}")));
      await session.step(64, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(65, "Then checkbox of \"bddmanaged\" membership row should be checked", () => shouldBe(page, el("checkbox of \"bddmanaged\" membership row"), "checked"));
      await session.step(66, "When user clicks on CANCEL button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(67, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
      await session.step(69, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A group added as a member is nested (Groups-15)", async () => {
      await session.step(72, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(73, "And user types \"BDD-GM-Child-{time}\" into membership search", () => typeInto(page, session.text("BDD-GM-Child-{time}"), el("membership search")));
      await session.step(74, "And user clicks on add button of \"BDD-GM-Child-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-GM-Child-{time}\" membership candidate"))));
      await session.step(75, "Then \"BDD-GM-Child-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Child-{time}\" membership row")), "visible"));
      await session.step(76, "When user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(77, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(78, "And \"BDD-GM-Child-{time}\" should be a member of \"BDD-GM-Group-{time}\" on the server", () => memberOnServer(page, session.text("BDD-GM-Child-{time}"), session.text("BDD-GM-Group-{time}")));
      await session.step(79, "And \"Members\" section in context panel should contain text \"BDD-GM-Child-{time}\"", () => shouldContainText(page, el("\"Members\" section in context panel"), session.text("BDD-GM-Child-{time}")));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
      await session.step(81, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Removing a member leaves the others (Groups-12)", async () => {
      await session.step(84, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(85, "Then \"bddmanaged\" membership row should be visible", () => shouldBe(page, el("\"bddmanaged\" membership row"), "visible"));
      await session.step(86, "When user clicks on remove button of \"bddmanaged\" membership row", () => clickOn(page, el("remove button of \"bddmanaged\" membership row")));
      await session.step(87, "Then \"bddmanaged\" membership row should be absent", () => shouldBe(page, el("\"bddmanaged\" membership row"), "absent"));
      await session.step(88, "When user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(89, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(90, "And \"bddmanaged\" should not be a member of \"BDD-GM-Group-{time}\" on the server", () => notMemberOnServer(page, "bddmanaged", session.text("BDD-GM-Group-{time}")));
      await session.step(91, "And \"BDD-GM-Child-{time}\" should be a member of \"BDD-GM-Group-{time}\" on the server", () => memberOnServer(page, session.text("BDD-GM-Child-{time}"), session.text("BDD-GM-Group-{time}")));
      await session.step(92, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(93, "Then no errors should have been logged", () => noErrors(page));
      await session.step(94, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A user's personal group is not in the Groups view (Groups-18)", async () => {
      await session.step(97, "When user types \"bddmanaged\" into gallery search", () => typeInto(page, "bddmanaged", el("gallery search")));
      await session.step(98, "Then gallery counter should have text \"0\"", () => shouldHaveText(page, el("gallery counter"), "0"));
      await session.step(99, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(100, "Then no errors should have been logged", () => noErrors(page));
      await session.step(101, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A group added as a member of its own member is refused as a circular reference", async () => {
      await session.step(105, "Given a group named \"BDD-GM-Outer-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Outer-{time}")));
      await session.step(106, "And a group named \"BDD-GM-Inner-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Inner-{time}")));
      await session.step(107, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(108, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(109, "And user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(110, "And user types \"BDD-GM-Outer-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Outer-{time}"), el("gallery search")));
      await session.step(111, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(112, "When user clicks on \"BDD-GM-Outer-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Outer-{time}\" link in gallery"))));
      await session.step(113, "Then the context panel should show \"BDD-GM-Outer-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Outer-{time}")));
      await session.step(114, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(115, "Then \"BDD-GM-Outer-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Outer-{time} members\" dialog")), "visible"));
      await session.step(116, "When user types \"BDD-GM-Inner-{time}\" into membership search", () => typeInto(page, session.text("BDD-GM-Inner-{time}"), el("membership search")));
      await session.step(117, "And user clicks on add button of \"BDD-GM-Inner-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-GM-Inner-{time}\" membership candidate"))));
      await session.step(118, "Then \"BDD-GM-Inner-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Inner-{time}\" membership row")), "visible"));
      await session.step(119, "When user clicks on SAVE button in \"BDD-GM-Outer-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Outer-{time} members\" dialog"))));
      await session.step(120, "Then the \"BDD-GM-Outer-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Outer-{time} members")));
      await session.step(121, "And \"BDD-GM-Inner-{time}\" should be a member of \"BDD-GM-Outer-{time}\" on the server", () => memberOnServer(page, session.text("BDD-GM-Inner-{time}"), session.text("BDD-GM-Outer-{time}")));
      await session.step(122, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(123, "And user types \"BDD-GM-Inner-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Inner-{time}"), el("gallery search")));
      await session.step(124, "And user clicks on \"BDD-GM-Inner-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Inner-{time}\" link in gallery"))));
      await session.step(125, "Then the context panel should show \"BDD-GM-Inner-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Inner-{time}")));
      await session.step(126, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(127, "Then \"BDD-GM-Inner-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Inner-{time} members\" dialog")), "visible"));
      await session.step(128, "When user types \"BDD-GM-Outer-{time}\" into membership search", () => typeInto(page, session.text("BDD-GM-Outer-{time}"), el("membership search")));
      await session.step(129, "And user clicks on add button of \"BDD-GM-Outer-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-GM-Outer-{time}\" membership candidate"))));
      await session.step(130, "Then \"BDD-GM-Outer-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Outer-{time}\" membership row")), "visible"));
      await session.step(131, "When user clicks on SAVE button in \"BDD-GM-Inner-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Inner-{time} members\" dialog"))));
      await session.step(132, "Then an error balloon containing \"Circular reference detected\" should have been shown", () => errorBalloonText(page, "Circular reference detected"));
      await session.step(133, "And \"BDD-GM-Outer-{time}\" should not be a member of \"BDD-GM-Inner-{time}\" on the server", () => notMemberOnServer(page, session.text("BDD-GM-Outer-{time}"), session.text("BDD-GM-Inner-{time}")));
      await session.step(134, "And \"BDD-GM-Inner-{time}\" should be a member of \"BDD-GM-Outer-{time}\" on the server", () => memberOnServer(page, session.text("BDD-GM-Inner-{time}"), session.text("BDD-GM-Outer-{time}")));
      await session.step(135, "When user clicks on CANCEL button in \"BDD-GM-Inner-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Inner-{time} members\" dialog"))));
      await session.step(136, "Then the \"BDD-GM-Inner-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Inner-{time} members")));
      await session.step(137, "When user clears gallery search", () => clearField(page, el("gallery search")));
    });
    await run.scenario("CANCEL in MANAGE leaves the members as they were", async () => {
      await session.step(141, "Given a group named \"BDD-GM-Edge-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Edge-{time}")));
      await session.step(142, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(143, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(144, "And user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(145, "And user types \"BDD-GM-Edge-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Edge-{time}"), el("gallery search")));
      await session.step(146, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(147, "When user clicks on \"BDD-GM-Edge-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Edge-{time}\" link in gallery"))));
      await session.step(148, "Then the context panel should show \"BDD-GM-Edge-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Edge-{time}")));
      await session.step(152, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(153, "Then \"BDD-GM-Edge-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Edge-{time} members\" dialog")), "visible"));
      await session.step(154, "When user types \"bddmanaged\" into membership search", () => typeInto(page, "bddmanaged", el("membership search")));
      await session.step(155, "And user clicks on add button of \"bddmanaged\" membership candidate", () => clickOn(page, el("add button of \"bddmanaged\" membership candidate")));
      await session.step(156, "Then \"bddmanaged\" membership row should be visible", () => shouldBe(page, el("\"bddmanaged\" membership row"), "visible"));
      await session.step(157, "When user clicks on CANCEL button in \"BDD-GM-Edge-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Edge-{time} members\" dialog"))));
      await session.step(158, "Then the \"BDD-GM-Edge-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Edge-{time} members")));
      await session.step(159, "And \"bddmanaged\" should not be a member of \"BDD-GM-Edge-{time}\" on the server", () => notMemberOnServer(page, "bddmanaged", session.text("BDD-GM-Edge-{time}")));
      await session.step(160, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(161, "Then no errors should have been logged", () => noErrors(page));
      await session.step(162, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The only admin unticking their Admin box is refused and stays the admin [change=unticking their Admin box, gesture=unchecks checkbox]", async () => {
      await session.step(141, "Given a group named \"BDD-GM-Edge-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Edge-{time}")));
      await session.step(142, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(143, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(144, "And user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(145, "And user types \"BDD-GM-Edge-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Edge-{time}"), el("gallery search")));
      await session.step(146, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(147, "When user clicks on \"BDD-GM-Edge-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Edge-{time}\" link in gallery"))));
      await session.step(148, "Then the context panel should show \"BDD-GM-Edge-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Edge-{time}")));
      await session.step(166, "Then \"{login}\" should be an admin member of \"BDD-GM-Edge-{time}\" on the server", () => adminMemberOnServer(page, session.text("{login}"), session.text("BDD-GM-Edge-{time}")));
      await session.step(167, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(168, "Then \"BDD-GM-Edge-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Edge-{time} members\" dialog")), "visible"));
      await session.step(169, "When user unchecks checkbox of \"{login}\" membership row", () => uncheck(page, el(session.text("checkbox of \"{login}\" membership row"))));
      await session.step(170, "And user clicks on SAVE button in \"BDD-GM-Edge-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Edge-{time} members\" dialog"))));
      await session.step(171, "Then an error balloon containing \"should have at least one admin user\" should have been shown", () => errorBalloonText(page, "should have at least one admin user"));
      await session.step(172, "And \"BDD-GM-Edge-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Edge-{time} members\" dialog")), "visible"));
      await session.step(173, "And \"{login}\" should be an admin member of \"BDD-GM-Edge-{time}\" on the server", () => adminMemberOnServer(page, session.text("{login}"), session.text("BDD-GM-Edge-{time}")));
      await session.step(174, "When user clicks on CANCEL button in \"BDD-GM-Edge-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Edge-{time} members\" dialog"))));
      await session.step(175, "Then the \"BDD-GM-Edge-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Edge-{time} members")));
      await session.step(176, "When user clears gallery search", () => clearField(page, el("gallery search")));
    });
    await run.scenario("The only admin removing themselves is refused and stays the admin [change=removing themselves, gesture=clicks on remove button]", async () => {
      await session.step(141, "Given a group named \"BDD-GM-Edge-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Edge-{time}")));
      await session.step(142, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(143, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(144, "And user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(145, "And user types \"BDD-GM-Edge-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Edge-{time}"), el("gallery search")));
      await session.step(146, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(147, "When user clicks on \"BDD-GM-Edge-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Edge-{time}\" link in gallery"))));
      await session.step(148, "Then the context panel should show \"BDD-GM-Edge-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Edge-{time}")));
      await session.step(166, "Then \"{login}\" should be an admin member of \"BDD-GM-Edge-{time}\" on the server", () => adminMemberOnServer(page, session.text("{login}"), session.text("BDD-GM-Edge-{time}")));
      await session.step(167, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(168, "Then \"BDD-GM-Edge-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Edge-{time} members\" dialog")), "visible"));
      await session.step(169, "When user clicks on remove button of \"{login}\" membership row", () => clickOn(page, el(session.text("remove button of \"{login}\" membership row"))));
      await session.step(170, "And user clicks on SAVE button in \"BDD-GM-Edge-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Edge-{time} members\" dialog"))));
      await session.step(171, "Then an error balloon containing \"should have at least one admin user\" should have been shown", () => errorBalloonText(page, "should have at least one admin user"));
      await session.step(172, "And \"BDD-GM-Edge-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Edge-{time} members\" dialog")), "visible"));
      await session.step(173, "And \"{login}\" should be an admin member of \"BDD-GM-Edge-{time}\" on the server", () => adminMemberOnServer(page, session.text("{login}"), session.text("BDD-GM-Edge-{time}")));
      await session.step(174, "When user clicks on CANCEL button in \"BDD-GM-Edge-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Edge-{time} members\" dialog"))));
      await session.step(175, "Then the \"BDD-GM-Edge-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Edge-{time} members")));
      await session.step(176, "When user clears gallery search", () => clearField(page, el("gallery search")));
    });
    run.finish();
  });
});
