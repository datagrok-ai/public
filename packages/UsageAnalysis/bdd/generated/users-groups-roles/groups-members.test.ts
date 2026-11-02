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
  test("A group's members", {tag: ["@journey", "@serial", "@groups", "@realizes:views.groups"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And a user \"bddmanaged\" is on the server", () => userOnServer(page, "bddmanaged"));
    await session.step(28, "And a group named \"BDD-GM-Group-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Group-{time}")));
    await session.step(29, "And a group named \"BDD-GM-Child-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Child-{time}")));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(32, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(33, "And user clicks on \"Platform > Groups\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Groups\" tree node inside browse tree")));
    await run.scenario("MANAGE adds a user to the group (Groups-11)", async () => {
      await session.step(36, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(37, "And user types \"BDD-GM-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Group-{time}"), el("gallery search")));
      await session.step(38, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(39, "When user clicks on \"BDD-GM-Group-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Group-{time}\" link in gallery"))));
      await session.step(40, "Then the context panel should show \"BDD-GM-Group-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Group-{time}")));
      await session.step(41, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(42, "And user types \"BDD-GM-Child-{time}\" into membership search", () => typeInto(page, session.text("BDD-GM-Child-{time}"), el("membership search")));
      await session.step(43, "And user clicks on add button of \"BDD-GM-Child-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-GM-Child-{time}\" membership candidate"))));
      await session.step(44, "Then \"BDD-GM-Child-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Child-{time}\" membership row")), "visible"));
      await session.step(45, "When user clicks on CANCEL button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(46, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(47, "And \"BDD-GM-Child-{time}\" should not be a member of \"BDD-GM-Group-{time}\" on the server", () => notMemberOnServer(page, session.text("BDD-GM-Child-{time}"), session.text("BDD-GM-Group-{time}")));
      await session.step(48, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(49, "Then \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(50, "When user types \"bddmanaged\" into membership search", () => typeInto(page, "bddmanaged", el("membership search")));
      await session.step(51, "And user clicks on add button of \"bddmanaged\" membership candidate", () => clickOn(page, el("add button of \"bddmanaged\" membership candidate")));
      await session.step(52, "Then \"bddmanaged\" membership row should be visible", () => shouldBe(page, el("\"bddmanaged\" membership row"), "visible"));
      await session.step(53, "When user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(54, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(55, "And \"bddmanaged\" should be a plain member of \"BDD-GM-Group-{time}\" on the server", () => plainMemberOnServer(page, "bddmanaged", session.text("BDD-GM-Group-{time}")));
      await session.step(56, "And \"Members\" section in context panel should contain text \"bddmanaged\"", () => shouldContainText(page, el("\"Members\" section in context panel"), "bddmanaged"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
      await session.step(58, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The only admin unticking their Admin box is refused and stays the admin [change=unticking their Admin box, gesture=unchecks checkbox]", async () => {
      await session.step(61, "Then \"{login}\" should be an admin member of \"BDD-GM-Group-{time}\" on the server", () => adminMemberOnServer(page, session.text("{login}"), session.text("BDD-GM-Group-{time}")));
      await session.step(62, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(63, "Then \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(64, "When user unchecks checkbox of \"{login}\" membership row", () => uncheck(page, el(session.text("checkbox of \"{login}\" membership row"))));
      await session.step(65, "And user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(66, "Then an error balloon containing \"should have at least one admin user\" should have been shown", () => errorBalloonText(page, "should have at least one admin user"));
      await session.step(67, "And \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(68, "And \"{login}\" should be an admin member of \"BDD-GM-Group-{time}\" on the server", () => adminMemberOnServer(page, session.text("{login}"), session.text("BDD-GM-Group-{time}")));
      await session.step(69, "When user clicks on CANCEL button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(70, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
    });
    await run.scenario("The only admin removing themselves is refused and stays the admin [change=removing themselves, gesture=clicks on remove button]", async () => {
      await session.step(61, "Then \"{login}\" should be an admin member of \"BDD-GM-Group-{time}\" on the server", () => adminMemberOnServer(page, session.text("{login}"), session.text("BDD-GM-Group-{time}")));
      await session.step(62, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(63, "Then \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(64, "When user clicks on remove button of \"{login}\" membership row", () => clickOn(page, el(session.text("remove button of \"{login}\" membership row"))));
      await session.step(65, "And user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(66, "Then an error balloon containing \"should have at least one admin user\" should have been shown", () => errorBalloonText(page, "should have at least one admin user"));
      await session.step(67, "And \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(68, "And \"{login}\" should be an admin member of \"BDD-GM-Group-{time}\" on the server", () => adminMemberOnServer(page, session.text("{login}"), session.text("BDD-GM-Group-{time}")));
      await session.step(69, "When user clicks on CANCEL button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(70, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
    });
    await run.scenario("The Admin box makes the member an admin (Groups-13)", async () => {
      await session.step(78, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(79, "Then checkbox label of \"bddmanaged\" membership row should have text \"Admin\"", () => shouldHaveText(page, el("checkbox label of \"bddmanaged\" membership row"), "Admin"));
      await session.step(80, "And checkbox of \"bddmanaged\" membership row should be unchecked", () => shouldBe(page, el("checkbox of \"bddmanaged\" membership row"), "unchecked"));
      await session.step(81, "When user checks checkbox of \"bddmanaged\" membership row", () => check(page, el("checkbox of \"bddmanaged\" membership row")));
      await session.step(82, "And user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(83, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(84, "And \"bddmanaged\" should be an admin member of \"BDD-GM-Group-{time}\" on the server", () => adminMemberOnServer(page, "bddmanaged", session.text("BDD-GM-Group-{time}")));
      await session.step(85, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(86, "Then checkbox of \"bddmanaged\" membership row should be checked", () => shouldBe(page, el("checkbox of \"bddmanaged\" membership row"), "checked"));
      await session.step(87, "When user clicks on CANCEL button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(88, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
      await session.step(90, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A group added as a member is nested (Groups-15)", async () => {
      await session.step(93, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(94, "And user types \"BDD-GM-Child-{time}\" into membership search", () => typeInto(page, session.text("BDD-GM-Child-{time}"), el("membership search")));
      await session.step(95, "And user clicks on add button of \"BDD-GM-Child-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-GM-Child-{time}\" membership candidate"))));
      await session.step(96, "Then \"BDD-GM-Child-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Child-{time}\" membership row")), "visible"));
      await session.step(97, "When user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(98, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(99, "And \"BDD-GM-Child-{time}\" should be a member of \"BDD-GM-Group-{time}\" on the server", () => memberOnServer(page, session.text("BDD-GM-Child-{time}"), session.text("BDD-GM-Group-{time}")));
      await session.step(100, "And \"Members\" section in context panel should contain text \"BDD-GM-Child-{time}\"", () => shouldContainText(page, el("\"Members\" section in context panel"), session.text("BDD-GM-Child-{time}")));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
      await session.step(102, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(103, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(104, "And user types \"BDD-GM-Child-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Child-{time}"), el("gallery search")));
      await session.step(105, "And user clicks on \"BDD-GM-Child-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Child-{time}\" link in gallery"))));
      await session.step(106, "Then the context panel should show \"BDD-GM-Child-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Child-{time}")));
      await session.step(107, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(108, "Then \"BDD-GM-Child-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Child-{time} members\" dialog")), "visible"));
      await session.step(109, "When user types \"BDD-GM-Group-{time}\" into membership search", () => typeInto(page, session.text("BDD-GM-Group-{time}"), el("membership search")));
      await session.step(110, "And user clicks on add button of \"BDD-GM-Group-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-GM-Group-{time}\" membership candidate"))));
      await session.step(111, "Then \"BDD-GM-Group-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time}\" membership row")), "visible"));
      await session.step(112, "When user clicks on SAVE button in \"BDD-GM-Child-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Child-{time} members\" dialog"))));
      await session.step(113, "Then an error balloon containing \"Circular reference detected\" should have been shown", () => errorBalloonText(page, "Circular reference detected"));
      await session.step(114, "And \"BDD-GM-Group-{time}\" should not be a member of \"BDD-GM-Child-{time}\" on the server", () => notMemberOnServer(page, session.text("BDD-GM-Group-{time}"), session.text("BDD-GM-Child-{time}")));
      await session.step(115, "And \"BDD-GM-Child-{time}\" should be a member of \"BDD-GM-Group-{time}\" on the server", () => memberOnServer(page, session.text("BDD-GM-Child-{time}"), session.text("BDD-GM-Group-{time}")));
      await session.step(116, "When user clicks on CANCEL button in \"BDD-GM-Child-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Child-{time} members\" dialog"))));
      await session.step(117, "Then the \"BDD-GM-Child-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Child-{time} members")));
      await session.step(118, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(119, "And user types \"BDD-GM-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Group-{time}"), el("gallery search")));
      await session.step(120, "And user clicks on \"BDD-GM-Group-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Group-{time}\" link in gallery"))));
      await session.step(121, "Then the context panel should show \"BDD-GM-Group-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Group-{time}")));
    });
    await run.scenario("Removing a member leaves the others (Groups-12)", async () => {
      await session.step(124, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(125, "Then \"bddmanaged\" membership row should be visible", () => shouldBe(page, el("\"bddmanaged\" membership row"), "visible"));
      await session.step(126, "When user clicks on remove button of \"bddmanaged\" membership row", () => clickOn(page, el("remove button of \"bddmanaged\" membership row")));
      await session.step(127, "Then \"bddmanaged\" membership row should be absent", () => shouldBe(page, el("\"bddmanaged\" membership row"), "absent"));
      await session.step(128, "When user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(129, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(130, "And \"bddmanaged\" should not be a member of \"BDD-GM-Group-{time}\" on the server", () => notMemberOnServer(page, "bddmanaged", session.text("BDD-GM-Group-{time}")));
      await session.step(131, "And \"BDD-GM-Child-{time}\" should be a member of \"BDD-GM-Group-{time}\" on the server", () => memberOnServer(page, session.text("BDD-GM-Child-{time}"), session.text("BDD-GM-Group-{time}")));
      await session.step(132, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(133, "Then no errors should have been logged", () => noErrors(page));
      await session.step(134, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A user's personal group is not in the Groups view (Groups-18)", async () => {
      await session.step(137, "When user types \"bddmanaged\" into gallery search", () => typeInto(page, "bddmanaged", el("gallery search")));
      await session.step(138, "Then gallery counter should have text \"0\"", () => shouldHaveText(page, el("gallery counter"), "0"));
      await session.step(139, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(140, "Then no errors should have been logged", () => noErrors(page));
      await session.step(141, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
