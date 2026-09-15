/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/groups-members.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.groups]
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
import {check, clearField, clickOn, expand, shouldBe, shouldContainText, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {adminMemberOnServer, browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, galleryCountLower, groupOnServer, memberOnServer, newUserOnServer, notMemberOnServer, plainMemberOnServer, rememberGalleryCount} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A group's members", () => {
  const session = feature(test, "features/users-groups-roles/groups-members.feature", import.meta.url);
  test("A group's members", {tag: ["@journey", "@groups", "@realizes:views.groups"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And a new user \"opavlenko{time}g\" with email \"opavlenko+{time}g@datagrok.ai\" is on the server", () => newUserOnServer(page, session.text("opavlenko{time}g"), session.text("opavlenko+{time}g@datagrok.ai")));
    await session.step(15, "And a group named \"BDD-GM-Group-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Group-{time}")));
    await session.step(16, "And a group named \"BDD-GM-Child-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GM-Child-{time}")));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(19, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(20, "And user clicks on \"Platform > Groups\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Groups\" tree node inside browse tree")));
    await run.scenario("MANAGE adds a user to the group (Groups-11)", async () => {
      await session.step(23, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(24, "And user types \"BDD-GM-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GM-Group-{time}"), el("gallery search")));
      await session.step(25, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(26, "When user clicks on \"BDD-GM-Group-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GM-Group-{time}\" link in gallery"))));
      await session.step(27, "Then the context panel should show \"BDD-GM-Group-{time}\"", () => contextPanelShows(page, session.text("BDD-GM-Group-{time}")));
      await session.step(28, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(29, "Then \"BDD-GM-Group-{time} members\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Group-{time} members\" dialog")), "visible"));
      await session.step(30, "When user types \"opavlenko{time}g\" into membership search", () => typeInto(page, session.text("opavlenko{time}g"), el("membership search")));
      await session.step(31, "And user clicks on add button of \"opavlenko{time}g\" membership candidate", () => clickOn(page, el(session.text("add button of \"opavlenko{time}g\" membership candidate"))));
      await session.step(32, "Then \"opavlenko{time}g\" membership row should be visible", () => shouldBe(page, el(session.text("\"opavlenko{time}g\" membership row")), "visible"));
      await session.step(33, "When user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(34, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(35, "And \"opavlenko{time}g\" should be a plain member of \"BDD-GM-Group-{time}\" on the server", () => plainMemberOnServer(page, session.text("opavlenko{time}g"), session.text("BDD-GM-Group-{time}")));
      await session.step(36, "And \"Members\" section in context panel should contain text \"opavlenko{time}g\"", () => shouldContainText(page, el("\"Members\" section in context panel"), session.text("opavlenko{time}g")));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
      await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Admin box makes the member an admin (Groups-13)", async () => {
      await session.step(41, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(42, "Then checkbox label of \"opavlenko{time}g\" membership row should have text \"Admin\"", () => shouldHaveText(page, el(session.text("checkbox label of \"opavlenko{time}g\" membership row")), "Admin"));
      await session.step(43, "And checkbox of \"opavlenko{time}g\" membership row should be unchecked", () => shouldBe(page, el(session.text("checkbox of \"opavlenko{time}g\" membership row")), "unchecked"));
      await session.step(44, "When user checks checkbox of \"opavlenko{time}g\" membership row", () => check(page, el(session.text("checkbox of \"opavlenko{time}g\" membership row"))));
      await session.step(45, "And user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(46, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(47, "And \"opavlenko{time}g\" should be an admin member of \"BDD-GM-Group-{time}\" on the server", () => adminMemberOnServer(page, session.text("opavlenko{time}g"), session.text("BDD-GM-Group-{time}")));
      await session.step(48, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(49, "Then checkbox of \"opavlenko{time}g\" membership row should be checked", () => shouldBe(page, el(session.text("checkbox of \"opavlenko{time}g\" membership row")), "checked"));
      await session.step(50, "When user clicks on CANCEL button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(51, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
      await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A group added as a member is nested (Groups-15)", async () => {
      await session.step(56, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(57, "And user types \"BDD-GM-Child-{time}\" into membership search", () => typeInto(page, session.text("BDD-GM-Child-{time}"), el("membership search")));
      await session.step(58, "And user clicks on add button of \"BDD-GM-Child-{time}\" membership candidate", () => clickOn(page, el(session.text("add button of \"BDD-GM-Child-{time}\" membership candidate"))));
      await session.step(59, "Then \"BDD-GM-Child-{time}\" membership row should be visible", () => shouldBe(page, el(session.text("\"BDD-GM-Child-{time}\" membership row")), "visible"));
      await session.step(60, "When user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(61, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(62, "And \"BDD-GM-Child-{time}\" should be a member of \"BDD-GM-Group-{time}\" on the server", () => memberOnServer(page, session.text("BDD-GM-Child-{time}"), session.text("BDD-GM-Group-{time}")));
      await session.step(63, "And \"Members\" section in context panel should contain text \"BDD-GM-Child-{time}\"", () => shouldContainText(page, el("\"Members\" section in context panel"), session.text("BDD-GM-Child-{time}")));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
      await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Removing a member leaves the others (Groups-12)", async () => {
      await session.step(68, "When user clicks on MANAGE button in \"Members\" section in context panel", () => clickOn(page, el("MANAGE button in \"Members\" section in context panel")));
      await session.step(69, "Then \"opavlenko{time}g\" membership row should be visible", () => shouldBe(page, el(session.text("\"opavlenko{time}g\" membership row")), "visible"));
      await session.step(70, "When user clicks on remove button of \"opavlenko{time}g\" membership row", () => clickOn(page, el(session.text("remove button of \"opavlenko{time}g\" membership row"))));
      await session.step(71, "Then \"opavlenko{time}g\" membership row should be absent", () => shouldBe(page, el(session.text("\"opavlenko{time}g\" membership row")), "absent"));
      await session.step(72, "When user clicks on SAVE button in \"BDD-GM-Group-{time} members\" dialog", () => clickOn(page, el(session.text("SAVE button in \"BDD-GM-Group-{time} members\" dialog"))));
      await session.step(73, "Then the \"BDD-GM-Group-{time} members\" dialog should close", () => dialogCloses(page, session.text("BDD-GM-Group-{time} members")));
      await session.step(74, "And \"opavlenko{time}g\" should not be a member of \"BDD-GM-Group-{time}\" on the server", () => notMemberOnServer(page, session.text("opavlenko{time}g"), session.text("BDD-GM-Group-{time}")));
      await session.step(75, "And \"BDD-GM-Child-{time}\" should be a member of \"BDD-GM-Group-{time}\" on the server", () => memberOnServer(page, session.text("BDD-GM-Child-{time}"), session.text("BDD-GM-Group-{time}")));
      await session.step(76, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(77, "Then no errors should have been logged", () => noErrors(page));
      await session.step(78, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
