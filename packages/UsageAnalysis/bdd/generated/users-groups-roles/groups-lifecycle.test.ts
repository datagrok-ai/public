/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/groups-lifecycle.feature
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
import {clearField, clickOn, expand, followingShouldBe, shouldBe, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, dialogCloses, galleryCountHigher, galleryCountLower, groupsOnServer, noGroupOnServer, rememberGalleryCount, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A group from creation to deletion", () => {
  const session = feature(test, "features/users-groups-roles/groups-lifecycle.feature", import.meta.url);
  test("A group from creation to deletion", {tag: ["@journey", "@serial", "@groups", "@realizes:views.groups"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And no group named \"BDD-GL-Group-{time}, BDD-GL-Renamed-{time}, BDD-GL-Cancelled-{time}\" is on the server", () => noGroupOnServer(page, session.text("BDD-GL-Group-{time}, BDD-GL-Renamed-{time}, BDD-GL-Cancelled-{time}")));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(20, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(21, "And user clicks on \"Platform > Groups\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Groups\" tree node inside browse tree")));
    await run.scenario("The New Group dialog, cancelled (Groups-03)", async () => {
      await session.step(24, "Then the \"Groups\" view should be current", () => viewIsCurrent(page, "Groups"));
      await session.step(25, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(26, "And user clicks on \"New Group...\" button", () => clickOn(page, el("\"New Group...\" button")));
      await session.step(27, "Then \"Create New Group\" dialog should be visible", () => shouldBe(page, el("\"Create New Group\" dialog"), "visible"));
      await session.step(28, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Name input in \"Create New Group\" dialog"],["Description input in \"Create New Group\" dialog"]]));
      await session.step(31, "When user types \"BDD-GL-Cancelled-{time}\" into Name input in \"Create New Group\" dialog", () => typeInto(page, session.text("BDD-GL-Cancelled-{time}"), el("Name input in \"Create New Group\" dialog")));
      await session.step(32, "And user clicks on CANCEL button in \"Create New Group\" dialog", () => clickOn(page, el("CANCEL button in \"Create New Group\" dialog")));
      await session.step(33, "Then the \"Create New Group\" dialog should close", () => dialogCloses(page, "Create New Group"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
      await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The New Group dialog refuses an empty name (Groups-05)", async () => {
      await session.step(38, "When user clicks on \"New Group...\" button", () => clickOn(page, el("\"New Group...\" button")));
      await session.step(39, "Then Name input in \"Create New Group\" dialog should have value \"New Group\"", () => shouldHaveValue(page, el("Name input in \"Create New Group\" dialog"), "New Group"));
      await session.step(40, "And OK button in \"Create New Group\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create New Group\" dialog"), "enabled"));
      await session.step(41, "When user clears Name input in \"Create New Group\" dialog", () => clearField(page, el("Name input in \"Create New Group\" dialog")));
      await session.step(42, "Then Name input in \"Create New Group\" dialog should be invalid", () => shouldBe(page, el("Name input in \"Create New Group\" dialog"), "invalid"));
      await session.step(43, "And OK button in \"Create New Group\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create New Group\" dialog"), "disabled"));
      await session.step(44, "When user clicks on CANCEL button in \"Create New Group\" dialog", () => clickOn(page, el("CANCEL button in \"Create New Group\" dialog")));
      await session.step(45, "Then the \"Create New Group\" dialog should close", () => dialogCloses(page, "Create New Group"));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A group is created (Groups-04)", async () => {
      await session.step(50, "When user clicks on \"New Group...\" button", () => clickOn(page, el("\"New Group...\" button")));
      await session.step(51, "And user types \"BDD-GL-Group-{time}\" into Name input in \"Create New Group\" dialog", () => typeInto(page, session.text("BDD-GL-Group-{time}"), el("Name input in \"Create New Group\" dialog")));
      await session.step(52, "And user types \"created by the BDD suite\" into Description input in \"Create New Group\" dialog", () => typeInto(page, "created by the BDD suite", el("Description input in \"Create New Group\" dialog")));
      await session.step(53, "And user clicks on OK button in \"Create New Group\" dialog", () => clickOn(page, el("OK button in \"Create New Group\" dialog")));
      await session.step(54, "Then the \"Create New Group\" dialog should close", () => dialogCloses(page, "Create New Group"));
      await session.step(55, "And 1 group named \"BDD-GL-Group-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-GL-Group-{time}")));
      await session.step(56, "And 0 groups named \"BDD-GL-Cancelled-{time}\" should be on the server", () => groupsOnServer(page, 0, session.text("BDD-GL-Cancelled-{time}")));
      await session.step(57, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(58, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
      await session.step(59, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(60, "And user types \"BDD-GL-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GL-Group-{time}"), el("gallery search")));
      await session.step(61, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(62, "And \"BDD-GL-Group-{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-GL-Group-{time}\" link in gallery")), "visible"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
      await session.step(64, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Properties... renames the group and rewrites its description (Groups-09)", async () => {
      await session.step(67, "When user picks \"Properties...\" from the context menu of \"BDD-GL-Group-{time}\" link in gallery", () => pickFromContextMenu(page, "Properties...", el(session.text("\"BDD-GL-Group-{time}\" link in gallery"))));
      await session.step(68, "Then \"BDD-GL-Group-{time} Properties\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GL-Group-{time} Properties\" dialog")), "visible"));
      await session.step(69, "And Name input in \"BDD-GL-Group-{time} Properties\" dialog should have value \"BDD-GL-Group-{time}\"", () => shouldHaveValue(page, el(session.text("Name input in \"BDD-GL-Group-{time} Properties\" dialog")), session.text("BDD-GL-Group-{time}")));
      await session.step(70, "And Description input in \"BDD-GL-Group-{time} Properties\" dialog should have value \"created by the BDD suite\"", () => shouldHaveValue(page, el(session.text("Description input in \"BDD-GL-Group-{time} Properties\" dialog")), "created by the BDD suite"));
      await session.step(71, "When user types \"BDD-GL-Renamed-{time}\" into Name input in \"BDD-GL-Group-{time} Properties\" dialog", () => typeInto(page, session.text("BDD-GL-Renamed-{time}"), el(session.text("Name input in \"BDD-GL-Group-{time} Properties\" dialog"))));
      await session.step(72, "And user types \"renamed by the BDD suite\" into Description input in \"BDD-GL-Group-{time} Properties\" dialog", () => typeInto(page, "renamed by the BDD suite", el(session.text("Description input in \"BDD-GL-Group-{time} Properties\" dialog"))));
      await session.step(73, "And user clicks on OK button in \"BDD-GL-Group-{time} Properties\" dialog", () => clickOn(page, el(session.text("OK button in \"BDD-GL-Group-{time} Properties\" dialog"))));
      await session.step(74, "Then the \"BDD-GL-Group-{time} Properties\" dialog should close", () => dialogCloses(page, session.text("BDD-GL-Group-{time} Properties")));
      await session.step(75, "And 1 group named \"BDD-GL-Renamed-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-GL-Renamed-{time}")));
      await session.step(76, "And 0 groups named \"BDD-GL-Group-{time}\" should be on the server", () => groupsOnServer(page, 0, session.text("BDD-GL-Group-{time}")));
      await session.step(77, "When user types \"BDD-GL-Renamed-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GL-Renamed-{time}"), el("gallery search")));
      await session.step(78, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(79, "When user picks \"Properties...\" from the context menu of \"BDD-GL-Renamed-{time}\" link in gallery", () => pickFromContextMenu(page, "Properties...", el(session.text("\"BDD-GL-Renamed-{time}\" link in gallery"))));
      await session.step(80, "Then \"BDD-GL-Renamed-{time} Properties\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GL-Renamed-{time} Properties\" dialog")), "visible"));
      await session.step(81, "And Description input in \"BDD-GL-Renamed-{time} Properties\" dialog should have value \"renamed by the BDD suite\"", () => shouldHaveValue(page, el(session.text("Description input in \"BDD-GL-Renamed-{time} Properties\" dialog")), "renamed by the BDD suite"));
      await session.step(82, "When user clicks on CANCEL button in \"BDD-GL-Renamed-{time} Properties\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GL-Renamed-{time} Properties\" dialog"))));
      await session.step(83, "Then the \"BDD-GL-Renamed-{time} Properties\" dialog should close", () => dialogCloses(page, session.text("BDD-GL-Renamed-{time} Properties")));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
      await session.step(85, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Delete removes the group after confirmation (Groups-14)", async () => {
      await session.step(88, "When user picks \"Delete\" from the context menu of \"BDD-GL-Renamed-{time}\" link in gallery", () => pickFromContextMenu(page, "Delete", el(session.text("\"BDD-GL-Renamed-{time}\" link in gallery"))));
      await session.step(89, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(90, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(91, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(92, "And 0 groups named \"BDD-GL-Renamed-{time}\" should be on the server", () => groupsOnServer(page, 0, session.text("BDD-GL-Renamed-{time}")));
      await session.step(93, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(94, "And user types \"BDD-GL-Renamed-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GL-Renamed-{time}"), el("gallery search")));
      await session.step(95, "Then \"BDD-GL-Renamed-{time}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BDD-GL-Renamed-{time}\" link in gallery")), "absent"));
      await session.step(96, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(97, "Then no errors should have been logged", () => noErrors(page));
      await session.step(98, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
