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
import {adminMemberOnServer, browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, galleryCountLower, galleryCountNotLower, groupOnServer, memberOnServer, notMemberOnServer, plainMemberOnServer, rememberGalleryCount, userOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A group's members", () => {
  const session = feature(test, "features/users-groups-roles/groups-members.feature", import.meta.url);
  test("A group's members", {tag: ["@journey", "@serial", "@groups", "@realizes:views.groups", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And a user \"bddmanaged\" is on the server", () => userOnServer(page, "bddmanaged"));
    await session.step(29, "And a group named \"BDD-GM-Group-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Group-{time}")));
    await session.step(30, "And a group named \"BDD-GM-Child-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Child-{time}")));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(32, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(33, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(34, "And user clicks on \"Platform > Groups\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Groups\" tree node inside browse tree")));
    await run.scenario("MANAGE adds a user to the group (Groups-11)", async () => {
      await session.step(37, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(38, "And user types \"BDD-GM-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Group-{time}"), el("gallery search")));
      await session.step(39, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(40, "When user clicks on \"BDD-GM-Group-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Group-{time}\" link in gallery"))));
      await session.step(41, "Then the context panel should show \"BDD-GM-Group-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Group-{time}")));
      await session.step(42, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(43, "Then \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(44, "When user types \"bddmanaged\" into membership search", () => typeInto(page, "bddmanaged", el("membership search")));
      await session.step(45, "And user clicks on add button of \"bddmanaged\" membership candidate", () => clickOn(page, el("add button of \"bddmanaged\" membership candidate")));
      await session.step(46, "Then \"bddmanaged\" membership row should be visible", () => shouldBe(page, el("\"bddmanaged\" membership row"), "visible"));
      await session.step(47, "When user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(48, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(49, "And \"bddmanaged\" should be a plain member of \"BDD-GM-Group-{time}\" on the server", () => plainMemberOnServer(page, "bddmanaged", session.text("BDD-GM-Group-{time}")));
      await session.step(50, "And \"Members\" section in context panel should contain text \"bddmanaged\"", () => shouldContainText(page, el("\"Members\" section in context panel"), "bddmanaged"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
      await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Admin box makes the member an admin (Groups-13)", async () => {
      await session.step(55, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(56, "Then checkbox label of \"bddmanaged\" membership row should have text \"Admin\"", () => shouldHaveText(page, el("checkbox label of \"bddmanaged\" membership row"), "Admin"));
      await session.step(57, "And checkbox of \"bddmanaged\" membership row should be unchecked", () => shouldBe(page, el("checkbox of \"bddmanaged\" membership row"), "unchecked"));
      await session.step(58, "When user checks checkbox of \"bddmanaged\" membership row", () => check(page, el("checkbox of \"bddmanaged\" membership row")));
      await session.step(59, "And user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(60, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(61, "And \"bddmanaged\" should be an admin member of \"BDD-GM-Group-{time}\" on the server", () => adminMemberOnServer(page, "bddmanaged", session.text("BDD-GM-Group-{time}")));
      await session.step(62, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(63, "Then checkbox of \"bddmanaged\" membership row should be checked", () => shouldBe(page, el("checkbox of \"bddmanaged\" membership row"), "checked"));
      await session.step(64, "When user clicks on CANCEL button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(65, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A group added as a member is nested (Groups-15)", async () => {
      await session.step(70, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(71, "And user types \"BDD-GM-Child-{time}\" into membership search", () => typeInto(page, session.text("BDD-GM-Child-{time}"), el("membership search")));
      await session.step(72, "And user clicks on add button of \"BDD-GM-Child-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-GM-Child-{time}\" membership candidate"))));
      await session.step(73, "Then \"BDD-GM-Child-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Child-{time}\" membership row")), "visible"));
      await session.step(74, "When user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(75, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(76, "And \"BDD-GM-Child-{time}\" should be a member of \"BDD-GM-Group-{time}\" on the server", () => memberOnServer(page, session.text("BDD-GM-Child-{time}"), session.text("BDD-GM-Group-{time}")));
      await session.step(77, "And \"Members\" section in context panel should contain text \"BDD-GM-Child-{time}\"", () => shouldContainText(page, el("\"Members\" section in context panel"), session.text("BDD-GM-Child-{time}")));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
      await session.step(79, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Removing a member leaves the others (Groups-12)", async () => {
      await session.step(82, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(83, "Then \"bddmanaged\" membership row should be visible", () => shouldBe(page, el("\"bddmanaged\" membership row"), "visible"));
      await session.step(84, "When user clicks on remove button of \"bddmanaged\" membership row", () => clickOn(page, el("remove button of \"bddmanaged\" membership row")));
      await session.step(85, "Then \"bddmanaged\" membership row should be absent", () => shouldBe(page, el("\"bddmanaged\" membership row"), "absent"));
      await session.step(86, "When user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(87, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(88, "And \"bddmanaged\" should not be a member of \"BDD-GM-Group-{time}\" on the server", () => notMemberOnServer(page, "bddmanaged", session.text("BDD-GM-Group-{time}")));
      await session.step(89, "And \"BDD-GM-Child-{time}\" should be a member of \"BDD-GM-Group-{time}\" on the server", () => memberOnServer(page, session.text("BDD-GM-Child-{time}"), session.text("BDD-GM-Group-{time}")));
      await session.step(90, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(91, "Then no errors should have been logged", () => noErrors(page));
      await session.step(92, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A user's personal group is not in the Groups view (Groups-18)", async () => {
      await session.step(95, "When user types \"bddmanaged\" into gallery search", () => typeInto(page, "bddmanaged", el("gallery search")));
      await session.step(96, "Then gallery counter should have text \"0\"", () => shouldHaveText(page, el("gallery counter"), "0"));
      await session.step(97, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(98, "Then no errors should have been logged", () => noErrors(page));
      await session.step(99, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("CANCEL in MANAGE leaves the members as they were", async () => {
      await session.step(103, "Then the gallery counter should not be lower than remembered", () => galleryCountNotLower(page));
      await session.step(104, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(105, "And user types \"BDD-GM-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Group-{time}"), el("gallery search")));
      await session.step(106, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(107, "When user clicks on \"BDD-GM-Group-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Group-{time}\" link in gallery"))));
      await session.step(108, "Then the context panel should show \"BDD-GM-Group-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Group-{time}")));
      await session.step(109, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(110, "Then \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(111, "When user types \"bddmanaged\" into membership search", () => typeInto(page, "bddmanaged", el("membership search")));
      await session.step(112, "And user clicks on add button of \"bddmanaged\" membership candidate", () => clickOn(page, el("add button of \"bddmanaged\" membership candidate")));
      await session.step(113, "Then \"bddmanaged\" membership row should be visible", () => shouldBe(page, el("\"bddmanaged\" membership row"), "visible"));
      await session.step(114, "When user clicks on CANCEL button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(115, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(116, "And \"bddmanaged\" should not be a member of \"BDD-GM-Group-{time}\" on the server", () => notMemberOnServer(page, "bddmanaged", session.text("BDD-GM-Group-{time}")));
      await session.step(117, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(118, "Then no errors should have been logged", () => noErrors(page));
      await session.step(119, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The only admin unticking their Admin box is refused and stays the admin [change=unticking their Admin box, gesture=unchecks checkbox]", async () => {
      await session.step(123, "Then the gallery counter should not be lower than remembered", () => galleryCountNotLower(page));
      await session.step(124, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(125, "And user types \"BDD-GM-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Group-{time}"), el("gallery search")));
      await session.step(126, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(127, "When user clicks on \"BDD-GM-Group-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Group-{time}\" link in gallery"))));
      await session.step(128, "Then the context panel should show \"BDD-GM-Group-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Group-{time}")));
      await session.step(129, "And \"{login}\" should be an admin member of \"BDD-GM-Group-{time}\" on the server", () => adminMemberOnServer(page, session.text("{login}"), session.text("BDD-GM-Group-{time}")));
      await session.step(130, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(131, "Then \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(132, "When user unchecks checkbox of \"{login}\" membership row", () => uncheck(page, el(session.text("checkbox of \"{login}\" membership row"))));
      await session.step(133, "And user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(134, "Then an error balloon containing \"should have at least one admin user\" should have been shown", () => errorBalloonText(page, "should have at least one admin user"));
      await session.step(135, "And \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(136, "And \"{login}\" should be an admin member of \"BDD-GM-Group-{time}\" on the server", () => adminMemberOnServer(page, session.text("{login}"), session.text("BDD-GM-Group-{time}")));
      await session.step(137, "When user clicks on CANCEL button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(138, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(139, "When user clears gallery search", () => clearField(page, el("gallery search")));
    });
    await run.scenario("The only admin removing themselves is refused and stays the admin [change=removing themselves, gesture=clicks on remove button]", async () => {
      await session.step(123, "Then the gallery counter should not be lower than remembered", () => galleryCountNotLower(page));
      await session.step(124, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(125, "And user types \"BDD-GM-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Group-{time}"), el("gallery search")));
      await session.step(126, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(127, "When user clicks on \"BDD-GM-Group-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Group-{time}\" link in gallery"))));
      await session.step(128, "Then the context panel should show \"BDD-GM-Group-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Group-{time}")));
      await session.step(129, "And \"{login}\" should be an admin member of \"BDD-GM-Group-{time}\" on the server", () => adminMemberOnServer(page, session.text("{login}"), session.text("BDD-GM-Group-{time}")));
      await session.step(130, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(131, "Then \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(132, "When user clicks on remove button of \"{login}\" membership row", () => clickOn(page, el(session.text("remove button of \"{login}\" membership row"))));
      await session.step(133, "And user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(134, "Then an error balloon containing \"should have at least one admin user\" should have been shown", () => errorBalloonText(page, "should have at least one admin user"));
      await session.step(135, "And \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(136, "And \"{login}\" should be an admin member of \"BDD-GM-Group-{time}\" on the server", () => adminMemberOnServer(page, session.text("{login}"), session.text("BDD-GM-Group-{time}")));
      await session.step(137, "When user clicks on CANCEL button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(138, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(139, "When user clears gallery search", () => clearField(page, el("gallery search")));
    });
    await run.scenario("A group added as a member of its own member is refused as a circular reference", async () => {
      await session.step(148, "Then the gallery counter should not be lower than remembered", () => galleryCountNotLower(page));
      await session.step(149, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(150, "And user types \"BDD-GM-Child-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Child-{time}"), el("gallery search")));
      await session.step(151, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(152, "When user clicks on \"BDD-GM-Child-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Child-{time}\" link in gallery"))));
      await session.step(153, "Then the context panel should show \"BDD-GM-Child-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Child-{time}")));
      await session.step(154, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(155, "Then \"BDD-GM-Child-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Child-{time} members\" dialog")), "visible"));
      await session.step(156, "When user types \"BDD-GM-Group-{time}\" into membership search", () => typeInto(page, session.text("BDD-GM-Group-{time}"), el("membership search")));
      await session.step(157, "And user clicks on add button of \"BDD-GM-Group-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-GM-Group-{time}\" membership candidate"))));
      await session.step(158, "Then \"BDD-GM-Group-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time}\" membership row")), "visible"));
      await session.step(159, "When user clicks on SAVE button in \"BDD-GM-Child-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Child-{time} members\" dialog"))));
      await session.step(160, "Then an error balloon containing \"Circular reference detected\" should have been shown", () => errorBalloonText(page, "Circular reference detected"));
      await session.step(161, "And \"BDD-GM-Group-{time}\" should not be a member of \"BDD-GM-Child-{time}\" on the server", () => notMemberOnServer(page, session.text("BDD-GM-Group-{time}"), session.text("BDD-GM-Child-{time}")));
      await session.step(162, "And \"BDD-GM-Child-{time}\" should be a member of \"BDD-GM-Group-{time}\" on the server", () => memberOnServer(page, session.text("BDD-GM-Child-{time}"), session.text("BDD-GM-Group-{time}")));
      await session.step(163, "When user clicks on CANCEL button in \"BDD-GM-Child-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Child-{time} members\" dialog"))));
      await session.step(164, "Then the \"BDD-GM-Child-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Child-{time} members")));
      await session.step(165, "When user clears gallery search", () => clearField(page, el("gallery search")));
    });
    run.finish();
  });
});
