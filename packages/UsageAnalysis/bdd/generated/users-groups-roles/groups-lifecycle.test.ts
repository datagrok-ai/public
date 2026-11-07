/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/groups-lifecycle.feature
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
import {clearField, clickOn, expand, followingShouldBe, shouldBe, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, dialogCloses, galleryCountHigher, galleryCountLower, groupOnServer, groupsOnServer, noGroupOnServer, rememberGalleryCount, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A group from creation to deletion", () => {
  const session = feature(test, "features/users-groups-roles/groups-lifecycle.feature", import.meta.url);
  test("A group from creation to deletion", {tag: ["@journey", "@serial", "@groups", "@realizes:views.groups", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And no group named \"BDD-GL-Group-{time}, BDD-GL-Renamed-{time}, BDD-GL-Cancelled-{time}\" is on the server", () => noGroupOnServer(page, session.text("BDD-GL-Group-{time}, BDD-GL-Renamed-{time}, BDD-GL-Cancelled-{time}")));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(26, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(27, "And user clicks on \"Platform > Groups\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Groups\" tree node inside browse tree")));
    await run.scenario("The New Group dialog, cancelled (Groups-03)", async () => {
      await session.step(30, "Then the \"Groups\" view should be current", () => viewIsCurrent(page, "Groups"));
      await session.step(31, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(32, "And user clicks on \"New Group...\" button", () => clickOn(page, el("\"New Group...\" button")));
      await session.step(33, "Then \"Create New Group\" dialog should be visible", () => shouldBe(page, el("\"Create New Group\" dialog"), "visible"));
      await session.step(34, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Name input in \"Create New Group\" dialog"],["Description input in \"Create New Group\" dialog"]]), [["Name input in \"Create New Group\" dialog"],["Description input in \"Create New Group\" dialog"]]);
      await session.step(37, "When user types \"BDD-GL-Cancelled-{time}\" into Name input in \"Create New Group\" dialog", () => typeInto(page, session.text("BDD-GL-Cancelled-{time}"), el("Name input in \"Create New Group\" dialog")));
      await session.step(38, "And user clicks on CANCEL button in \"Create New Group\" dialog", () => clickOn(page, el("CANCEL button in \"Create New Group\" dialog")));
      await session.step(39, "Then the \"Create New Group\" dialog should close", () => dialogCloses(page, "Create New Group"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
      await session.step(41, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The New Group dialog refuses an empty name (Groups-05)", async () => {
      await session.step(44, "When user clicks on \"New Group...\" button", () => clickOn(page, el("\"New Group...\" button")));
      await session.step(45, "Then Name input in \"Create New Group\" dialog should have value \"New Group\"", () => shouldHaveValue(page, el("Name input in \"Create New Group\" dialog"), "New Group"));
      await session.step(46, "And OK button in \"Create New Group\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create New Group\" dialog"), "enabled"));
      await session.step(47, "When user clears Name input in \"Create New Group\" dialog", () => clearField(page, el("Name input in \"Create New Group\" dialog")));
      await session.step(48, "Then Name input in \"Create New Group\" dialog should be invalid", () => shouldBe(page, el("Name input in \"Create New Group\" dialog"), "invalid"));
      await session.step(49, "And OK button in \"Create New Group\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create New Group\" dialog"), "disabled"));
      await session.step(50, "When user clicks on CANCEL button in \"Create New Group\" dialog", () => clickOn(page, el("CANCEL button in \"Create New Group\" dialog")));
      await session.step(51, "Then the \"Create New Group\" dialog should close", () => dialogCloses(page, "Create New Group"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
      await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A group is created (Groups-04)", async () => {
      await session.step(56, "When user clicks on \"New Group...\" button", () => clickOn(page, el("\"New Group...\" button")));
      await session.step(57, "And user types \"BDD-GL-Group-{time}\" into Name input in \"Create New Group\" dialog", () => typeInto(page, session.text("BDD-GL-Group-{time}"), el("Name input in \"Create New Group\" dialog")));
      await session.step(58, "And user types \"created by the BDD suite\" into Description input in \"Create New Group\" dialog", () => typeInto(page, "created by the BDD suite", el("Description input in \"Create New Group\" dialog")));
      await session.step(59, "And user clicks on OK button in \"Create New Group\" dialog", () => clickOn(page, el("OK button in \"Create New Group\" dialog")));
      await session.step(60, "Then the \"Create New Group\" dialog should close", () => dialogCloses(page, "Create New Group"));
      await session.step(61, "And 1 group named \"BDD-GL-Group-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-GL-Group-{time}")));
      await session.step(62, "And 0 groups named \"BDD-GL-Cancelled-{time}\" should be on the server", () => groupsOnServer(page, 0, session.text("BDD-GL-Cancelled-{time}")));
      await session.step(63, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(64, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
      await session.step(65, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(66, "And user types \"BDD-GL-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GL-Group-{time}"), el("gallery search")));
      await session.step(67, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(68, "And \"BDD-GL-Group-{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-GL-Group-{time}\" link in gallery")), "visible"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
      await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Properties... renames the group and rewrites its description (Groups-09)", async () => {
      await session.step(73, "When user picks \"Properties...\" from the context menu of \"BDD-GL-Group-{time}\" link in gallery", () => pickFromContextMenu(page, "Properties...", el(session.text("\"BDD-GL-Group-{time}\" link in gallery"))));
      await session.step(74, "Then \"BDD-GL-Group-{time} Properties\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GL-Group-{time} Properties\" dialog")), "visible"));
      await session.step(75, "And Name input in \"BDD-GL-Group-{time} Properties\" dialog should have value \"BDD-GL-Group-{time}\"", () => shouldHaveValue(page, el(session.text("Name input in \"BDD-GL-Group-{time} Properties\" dialog")), session.text("BDD-GL-Group-{time}")));
      await session.step(76, "And Description input in \"BDD-GL-Group-{time} Properties\" dialog should have value \"created by the BDD suite\"", () => shouldHaveValue(page, el(session.text("Description input in \"BDD-GL-Group-{time} Properties\" dialog")), "created by the BDD suite"));
      await session.step(77, "When user types \"BDD-GL-Renamed-{time}\" into Name input in \"BDD-GL-Group-{time} Properties\" dialog", () => typeInto(page, session.text("BDD-GL-Renamed-{time}"), el(session.text("Name input in \"BDD-GL-Group-{time} Properties\" dialog"))));
      await session.step(78, "And user types \"renamed by the BDD suite\" into Description input in \"BDD-GL-Group-{time} Properties\" dialog", () => typeInto(page, "renamed by the BDD suite", el(session.text("Description input in \"BDD-GL-Group-{time} Properties\" dialog"))));
      await session.step(79, "And user clicks on OK button in \"BDD-GL-Group-{time} Properties\" dialog", () => clickOn(page, el(session.text("OK button in \"BDD-GL-Group-{time} Properties\" dialog"))));
      await session.step(80, "Then the \"BDD-GL-Group-{time} Properties\" dialog should close", () => dialogCloses(page, session.text("BDD-GL-Group-{time} Properties")));
      await session.step(81, "And 1 group named \"BDD-GL-Renamed-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-GL-Renamed-{time}")));
      await session.step(82, "And 0 groups named \"BDD-GL-Group-{time}\" should be on the server", () => groupsOnServer(page, 0, session.text("BDD-GL-Group-{time}")));
      await session.step(83, "When user types \"BDD-GL-Renamed-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GL-Renamed-{time}"), el("gallery search")));
      await session.step(84, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(85, "When user picks \"Properties...\" from the context menu of \"BDD-GL-Renamed-{time}\" link in gallery", () => pickFromContextMenu(page, "Properties...", el(session.text("\"BDD-GL-Renamed-{time}\" link in gallery"))));
      await session.step(86, "Then \"BDD-GL-Renamed-{time} Properties\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-GL-Renamed-{time} Properties\" dialog")), "visible"));
      await session.step(87, "And Description input in \"BDD-GL-Renamed-{time} Properties\" dialog should have value \"renamed by the BDD suite\"", () => shouldHaveValue(page, el(session.text("Description input in \"BDD-GL-Renamed-{time} Properties\" dialog")), "renamed by the BDD suite"));
      await session.step(88, "When user clicks on CANCEL button in \"BDD-GL-Renamed-{time} Properties\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-GL-Renamed-{time} Properties\" dialog"))));
      await session.step(89, "Then the \"BDD-GL-Renamed-{time} Properties\" dialog should close", () => dialogCloses(page, session.text("BDD-GL-Renamed-{time} Properties")));
      await session.step(90, "And no errors should have been logged", () => noErrors(page));
      await session.step(91, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Delete removes the group after confirmation (Groups-14)", async () => {
      await session.step(94, "When user picks \"Delete\" from the context menu of \"BDD-GL-Renamed-{time}\" link in gallery", () => pickFromContextMenu(page, "Delete", el(session.text("\"BDD-GL-Renamed-{time}\" link in gallery"))));
      await session.step(95, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(96, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(97, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(98, "And 0 groups named \"BDD-GL-Renamed-{time}\" should be on the server", () => groupsOnServer(page, 0, session.text("BDD-GL-Renamed-{time}")));
      await session.step(99, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(100, "And user types \"BDD-GL-Renamed-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GL-Renamed-{time}"), el("gallery search")));
      await session.step(101, "Then \"BDD-GL-Renamed-{time}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BDD-GL-Renamed-{time}\" link in gallery")), "absent"));
      await session.step(102, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(103, "Then no errors should have been logged", () => noErrors(page));
      await session.step(104, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The New Group dialog refuses a name a group already has", async () => {
      await session.step(108, "Given a group named \"BDDGEApi{time}\" is on the server", () => groupOnServer(page, session.text("BDDGEApi{time}")));
      await session.step(109, "When user clicks on \"New Group...\" button", () => clickOn(page, el("\"New Group...\" button")));
      await session.step(110, "And user types \"BDDGEApi{time}\" into Name input in \"Create New Group\" dialog", () => typeInto(page, session.text("BDDGEApi{time}"), el("Name input in \"Create New Group\" dialog")));
      await session.step(111, "And user clicks on OK button in \"Create New Group\" dialog", () => clickOn(page, el("OK button in \"Create New Group\" dialog")));
      await session.step(112, "Then an error balloon containing \"A group or role with the same name already exists\" should have been shown", () => errorBalloonText(page, "A group or role with the same name already exists"));
      await session.step(113, "And \"Create New Group\" dialog should be visible", () => shouldBe(page, el("\"Create New Group\" dialog"), "visible"));
      await session.step(114, "And 1 group named \"BDDGEApi{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDDGEApi{time}")));
      await session.step(115, "When user clicks on CANCEL button in \"Create New Group\" dialog", () => clickOn(page, el("CANCEL button in \"Create New Group\" dialog")));
      await session.step(116, "Then the \"Create New Group\" dialog should close", () => dialogCloses(page, "Create New Group"));
    });
    await run.scenario("Properties... refuses a name of spaces, and CANCEL keeps the group's name", async () => {
      await session.step(120, "Given a group named \"BDDGEApi{time}\" is on the server", () => groupOnServer(page, session.text("BDDGEApi{time}")));
      await session.step(121, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(122, "And user types \"BDDGEApi{time}\" into gallery search", () => typeInto(page, session.text("BDDGEApi{time}"), el("gallery search")));
      await session.step(123, "Then \"BDDGEApi{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDGEApi{time}\" link in gallery")), "visible"));
      await session.step(124, "When user picks \"Properties...\" from the context menu of \"BDDGEApi{time}\" link in gallery", () => pickFromContextMenu(page, "Properties...", el(session.text("\"BDDGEApi{time}\" link in gallery"))));
      await session.step(125, "Then \"BDDGEApi{time} Properties\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDGEApi{time} Properties\" dialog")), "visible"));
      await session.step(126, "When user types \"   \" into Name input in \"BDDGEApi{time} Properties\" dialog", () => typeInto(page, "   ", el(session.text("Name input in \"BDDGEApi{time} Properties\" dialog"))));
      await session.step(127, "Then Name input in \"BDDGEApi{time} Properties\" dialog should be invalid", () => shouldBe(page, el(session.text("Name input in \"BDDGEApi{time} Properties\" dialog")), "invalid"));
      await session.step(128, "And OK button in \"BDDGEApi{time} Properties\" dialog should be disabled", () => shouldBe(page, el(session.text("OK button in \"BDDGEApi{time} Properties\" dialog")), "disabled"));
      await session.step(129, "When user clicks on CANCEL button in \"BDDGEApi{time} Properties\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDDGEApi{time} Properties\" dialog"))));
      await session.step(130, "Then the \"BDDGEApi{time} Properties\" dialog should close", () => dialogCloses(page, session.text("BDDGEApi{time} Properties")));
      await session.step(131, "And 1 group named \"BDDGEApi{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDDGEApi{time}")));
      await session.step(132, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(133, "Then no errors should have been logged", () => noErrors(page));
      await session.step(134, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
