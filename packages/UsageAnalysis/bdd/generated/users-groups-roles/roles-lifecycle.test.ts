/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/roles-lifecycle.feature
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
import {clearField, clickOn, expand, followingShouldBe, shouldBe, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, dialogCloses, galleryCountHigher, galleryCountLower, groupsOnServer, noGroupOnServer, rememberGalleryCount, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A role from creation to deletion", () => {
  const session = feature(test, "features/users-groups-roles/roles-lifecycle.feature", import.meta.url);
  test("A role from creation to deletion", {tag: ["@journey", "@serial", "@roles", "@realizes:views.roles"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And no role named \"BDD-RL-Role-{time}, BDD-RL-Renamed-{time}, BDD-RL-Cancelled-{time}\" is on the server", () => noGroupOnServer(page, session.text("BDD-RL-Role-{time}, BDD-RL-Renamed-{time}, BDD-RL-Cancelled-{time}")));
    await session.step(24, "And no group named \"BDD-RE-Group-{time}\" is on the server", () => noGroupOnServer(page, session.text("BDD-RE-Group-{time}")));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(27, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(28, "And user clicks on \"Platform > Roles\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Roles\" tree node inside browse tree")));
    await run.scenario("The New Role dialog, cancelled (Roles-03)", async () => {
      await session.step(31, "Then the \"Roles\" view should be current", () => viewIsCurrent(page, "Roles"));
      await session.step(32, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(33, "And user clicks on \"New Role...\" button", () => clickOn(page, el("\"New Role...\" button")));
      await session.step(34, "Then \"Create New Role\" dialog should be visible", () => shouldBe(page, el("\"Create New Role\" dialog"), "visible"));
      await session.step(35, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Name input in \"Create New Role\" dialog"],["Description input in \"Create New Role\" dialog"]]), [["Name input in \"Create New Role\" dialog"],["Description input in \"Create New Role\" dialog"]]);
      await session.step(38, "When user types \"BDD-RL-Cancelled-{time}\" into Name input in \"Create New Role\" dialog", () => typeInto(page, session.text("BDD-RL-Cancelled-{time}"), el("Name input in \"Create New Role\" dialog")));
      await session.step(39, "And user clicks on CANCEL button in \"Create New Role\" dialog", () => clickOn(page, el("CANCEL button in \"Create New Role\" dialog")));
      await session.step(40, "Then the \"Create New Role\" dialog should close", () => dialogCloses(page, "Create New Role"));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The New Role dialog refuses an empty name (Roles-05)", async () => {
      await session.step(45, "When user clicks on \"New Role...\" button", () => clickOn(page, el("\"New Role...\" button")));
      await session.step(46, "Then Name input in \"Create New Role\" dialog should have value \"New Role\"", () => shouldHaveValue(page, el("Name input in \"Create New Role\" dialog"), "New Role"));
      await session.step(47, "And OK button in \"Create New Role\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create New Role\" dialog"), "enabled"));
      await session.step(48, "When user types \"   \" into Name input in \"Create New Role\" dialog", () => typeInto(page, "   ", el("Name input in \"Create New Role\" dialog")));
      await session.step(49, "Then Name input in \"Create New Role\" dialog should be invalid", () => shouldBe(page, el("Name input in \"Create New Role\" dialog"), "invalid"));
      await session.step(50, "And OK button in \"Create New Role\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create New Role\" dialog"), "disabled"));
      await session.step(51, "When user types \"New Role\" into Name input in \"Create New Role\" dialog", () => typeInto(page, "New Role", el("Name input in \"Create New Role\" dialog")));
      await session.step(52, "Then OK button in \"Create New Role\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create New Role\" dialog"), "enabled"));
      await session.step(53, "When user clears Name input in \"Create New Role\" dialog", () => clearField(page, el("Name input in \"Create New Role\" dialog")));
      await session.step(54, "Then Name input in \"Create New Role\" dialog should be invalid", () => shouldBe(page, el("Name input in \"Create New Role\" dialog"), "invalid"));
      await session.step(55, "And OK button in \"Create New Role\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create New Role\" dialog"), "disabled"));
      await session.step(56, "When user clicks on CANCEL button in \"Create New Role\" dialog", () => clickOn(page, el("CANCEL button in \"Create New Role\" dialog")));
      await session.step(57, "Then the \"Create New Role\" dialog should close", () => dialogCloses(page, "Create New Role"));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
      await session.step(59, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(60, "When user clicks on \"Platform > Groups\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Groups\" tree node inside browse tree")));
      await session.step(61, "Then the \"Groups\" view should be current", () => viewIsCurrent(page, "Groups"));
      await session.step(62, "When user clicks on \"New Group...\" button", () => clickOn(page, el("\"New Group...\" button")));
      await session.step(63, "And user types \"BDD-RE-Group-{time}\" into Name input in \"Create New Group\" dialog", () => typeInto(page, session.text("BDD-RE-Group-{time}"), el("Name input in \"Create New Group\" dialog")));
      await session.step(64, "And user clicks on OK button in \"Create New Group\" dialog", () => clickOn(page, el("OK button in \"Create New Group\" dialog")));
      await session.step(65, "Then the \"Create New Group\" dialog should close", () => dialogCloses(page, "Create New Group"));
      await session.step(66, "And 1 group named \"BDD-RE-Group-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-RE-Group-{time}")));
      await session.step(67, "When user clicks on \"Platform > Roles\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Roles\" tree node inside browse tree")));
      await session.step(68, "Then the \"Roles\" view should be current", () => viewIsCurrent(page, "Roles"));
      await session.step(69, "When user clicks on \"New Role...\" button", () => clickOn(page, el("\"New Role...\" button")));
      await session.step(70, "And user types \"BDD-RE-Group-{time}\" into Name input in \"Create New Role\" dialog", () => typeInto(page, session.text("BDD-RE-Group-{time}"), el("Name input in \"Create New Role\" dialog")));
      await session.step(71, "And user clicks on OK button in \"Create New Role\" dialog", () => clickOn(page, el("OK button in \"Create New Role\" dialog")));
      await session.step(72, "Then an error balloon containing \"A group or role with the same name already exists\" should have been shown", () => errorBalloonText(page, "A group or role with the same name already exists"));
      await session.step(73, "And \"Create New Role\" dialog should be visible", () => shouldBe(page, el("\"Create New Role\" dialog"), "visible"));
      await session.step(74, "And 1 group named \"BDD-RE-Group-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-RE-Group-{time}")));
      await session.step(75, "When user clicks on CANCEL button in \"Create New Role\" dialog", () => clickOn(page, el("CANCEL button in \"Create New Role\" dialog")));
      await session.step(76, "Then the \"Create New Role\" dialog should close", () => dialogCloses(page, "Create New Role"));
    });
    await run.scenario("A role is created (Roles-04)", async () => {
      await session.step(79, "When user clicks on \"New Role...\" button", () => clickOn(page, el("\"New Role...\" button")));
      await session.step(80, "And user types \"BDD-RL-Role-{time}\" into Name input in \"Create New Role\" dialog", () => typeInto(page, session.text("BDD-RL-Role-{time}"), el("Name input in \"Create New Role\" dialog")));
      await session.step(81, "And user types \"created by the BDD suite\" into Description input in \"Create New Role\" dialog", () => typeInto(page, "created by the BDD suite", el("Description input in \"Create New Role\" dialog")));
      await session.step(82, "And user clicks on OK button in \"Create New Role\" dialog", () => clickOn(page, el("OK button in \"Create New Role\" dialog")));
      await session.step(83, "Then the \"Create New Role\" dialog should close", () => dialogCloses(page, "Create New Role"));
      await session.step(84, "And 1 role named \"BDD-RL-Role-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-RL-Role-{time}")));
      await session.step(85, "And 0 roles named \"BDD-RL-Cancelled-{time}\" should be on the server", () => groupsOnServer(page, 0, session.text("BDD-RL-Cancelled-{time}")));
      await session.step(86, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(87, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
      await session.step(88, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(89, "And user types \"BDD-RL-Role-{time}\" into gallery search", () => typeInto(page, session.text("BDD-RL-Role-{time}"), el("gallery search")));
      await session.step(90, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(91, "And \"BDD-RL-Role-{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-RL-Role-{time}\" link in gallery")), "visible"));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
      await session.step(93, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Properties... renames the role and rewrites its description (Roles-09)", async () => {
      await session.step(96, "When user picks \"Properties...\" from the context menu of \"BDD-RL-Role-{time}\" link in gallery", () => pickFromContextMenu(page, "Properties...", el(session.text("\"BDD-RL-Role-{time}\" link in gallery"))));
      await session.step(97, "Then \"BDD-RL-Role-{time} Properties\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-RL-Role-{time} Properties\" dialog")), "visible"));
      await session.step(98, "And Name input in \"BDD-RL-Role-{time} Properties\" dialog should have value \"BDD-RL-Role-{time}\"", () => shouldHaveValue(page, el(session.text("Name input in \"BDD-RL-Role-{time} Properties\" dialog")), session.text("BDD-RL-Role-{time}")));
      await session.step(99, "And Description input in \"BDD-RL-Role-{time} Properties\" dialog should have value \"created by the BDD suite\"", () => shouldHaveValue(page, el(session.text("Description input in \"BDD-RL-Role-{time} Properties\" dialog")), "created by the BDD suite"));
      await session.step(100, "When user types \"BDD-RL-Renamed-{time}\" into Name input in \"BDD-RL-Role-{time} Properties\" dialog", () => typeInto(page, session.text("BDD-RL-Renamed-{time}"), el(session.text("Name input in \"BDD-RL-Role-{time} Properties\" dialog"))));
      await session.step(101, "And user types \"renamed by the BDD suite\" into Description input in \"BDD-RL-Role-{time} Properties\" dialog", () => typeInto(page, "renamed by the BDD suite", el(session.text("Description input in \"BDD-RL-Role-{time} Properties\" dialog"))));
      await session.step(102, "And user clicks on OK button in \"BDD-RL-Role-{time} Properties\" dialog", () => clickOn(page, el(session.text("OK button in \"BDD-RL-Role-{time} Properties\" dialog"))));
      await session.step(103, "Then the \"BDD-RL-Role-{time} Properties\" dialog should close", () => dialogCloses(page, session.text("BDD-RL-Role-{time} Properties")));
      await session.step(104, "And 1 role named \"BDD-RL-Renamed-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-RL-Renamed-{time}")));
      await session.step(105, "And 0 roles named \"BDD-RL-Role-{time}\" should be on the server", () => groupsOnServer(page, 0, session.text("BDD-RL-Role-{time}")));
      await session.step(106, "When user types \"BDD-RL-Renamed-{time}\" into gallery search", () => typeInto(page, session.text("BDD-RL-Renamed-{time}"), el("gallery search")));
      await session.step(107, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(108, "When user picks \"Properties...\" from the context menu of \"BDD-RL-Renamed-{time}\" link in gallery", () => pickFromContextMenu(page, "Properties...", el(session.text("\"BDD-RL-Renamed-{time}\" link in gallery"))));
      await session.step(109, "Then \"BDD-RL-Renamed-{time} Properties\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDD-RL-Renamed-{time} Properties\" dialog")), "visible"));
      await session.step(110, "And Description input in \"BDD-RL-Renamed-{time} Properties\" dialog should have value \"renamed by the BDD suite\"", () => shouldHaveValue(page, el(session.text("Description input in \"BDD-RL-Renamed-{time} Properties\" dialog")), "renamed by the BDD suite"));
      await session.step(111, "When user clicks on CANCEL button in \"BDD-RL-Renamed-{time} Properties\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDD-RL-Renamed-{time} Properties\" dialog"))));
      await session.step(112, "Then the \"BDD-RL-Renamed-{time} Properties\" dialog should close", () => dialogCloses(page, session.text("BDD-RL-Renamed-{time} Properties")));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
      await session.step(114, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Delete removes the role after confirmation (Roles-15)", async () => {
      await session.step(117, "When user picks \"Delete\" from the context menu of \"BDD-RL-Renamed-{time}\" link in gallery", () => pickFromContextMenu(page, "Delete", el(session.text("\"BDD-RL-Renamed-{time}\" link in gallery"))));
      await session.step(118, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(119, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(120, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(121, "And 0 roles named \"BDD-RL-Renamed-{time}\" should be on the server", () => groupsOnServer(page, 0, session.text("BDD-RL-Renamed-{time}")));
      await session.step(122, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(123, "And user types \"BDD-RL-Renamed-{time}\" into gallery search", () => typeInto(page, session.text("BDD-RL-Renamed-{time}"), el("gallery search")));
      await session.step(124, "Then \"BDD-RL-Renamed-{time}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BDD-RL-Renamed-{time}\" link in gallery")), "absent"));
      await session.step(125, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(126, "Then no errors should have been logged", () => noErrors(page));
      await session.step(127, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
