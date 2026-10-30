/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/files-lifecycle.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
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
import {clickOn, enterInto, isExpanded, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses, refreshBrowse} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fixtureInUsersFiles} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Files: a folder and a file created, renamed and deleted from the Browse tree", () => {
  const session = feature(test, "features/browse/files-lifecycle.feature", import.meta.url);
  test("Files: a folder and a file created, renamed and deleted from the Browse tree", {tag: ["@journey", "@browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the \"fixtures/cars-small.csv\" file of the project is in the user's files as \"bdd-files-{time}/seed.csv\"", () => fixtureInUsersFiles(page, "fixtures/cars-small.csv", session.text("bdd-files-{time}/seed.csv")));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(26, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(27, "And Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
    await session.step(28, "And Files---My-files---bdd-files-{time} tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("Files---My-files---bdd-files-{time} tree node inside browse tree"))));
    await run.scenario("Create folder... asks for a name and the tree shows the new folder", async () => {
      await session.step(31, "Then Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree")), "visible"));
      await session.step(32, "When user picks \"Create folder...\" from the context menu of Files---My-files---bdd-files-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Create folder...", el(session.text("Files---My-files---bdd-files-{time} tree node inside browse tree"))));
      await session.step(33, "Then \"Create folder\" dialog should be visible", () => shouldBe(page, el("\"Create folder\" dialog"), "visible"));
      await session.step(34, "When user enters \"sub-a\" into Name input in \"Create folder\" dialog", () => enterInto(page, "sub-a", el("Name input in \"Create folder\" dialog")));
      await session.step(35, "And user clicks on OK button in \"Create folder\" dialog", () => clickOn(page, el("OK button in \"Create folder\" dialog")));
      await session.step(36, "Then the \"Create folder\" dialog should close", () => dialogCloses(page, "Create folder"));
      await session.step(37, "And Files---My-files---bdd-files-{time}---sub-a tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---sub-a tree node inside browse tree")), "visible"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Rename... shows the folder's name and the tree shows the new one", async () => {
      await session.step(41, "When user picks \"Rename...\" from the context menu of Files---My-files---bdd-files-{time}---sub-a tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Files---My-files---bdd-files-{time}---sub-a tree node inside browse tree"))));
      await session.step(42, "Then \"Folder name\" input in Rename dialog should have value \"sub-a\"", () => shouldHaveValue(page, el("\"Folder name\" input in Rename dialog"), "sub-a"));
      await session.step(43, "When user enters \"sub-b\" into \"Folder name\" input in Rename dialog", () => enterInto(page, "sub-b", el("\"Folder name\" input in Rename dialog")));
      await session.step(44, "And user clicks on OK button in Rename dialog", () => clickOn(page, el("OK button in Rename dialog")));
      await session.step(45, "Then the \"Rename\" dialog should close", () => dialogCloses(page, "Rename"));
      await session.step(46, "And Files---My-files---bdd-files-{time}---sub-b tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---sub-b tree node inside browse tree")), "visible"));
      await session.step(47, "And Files---My-files---bdd-files-{time}---sub-a tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---sub-a tree node inside browse tree")), "absent"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Delete... asks for confirmation and the folder leaves the tree", async () => {
      await session.step(51, "When user picks \"Delete...\" from the context menu of Files---My-files---bdd-files-{time}---sub-b tree node inside browse tree", () => pickFromContextMenu(page, "Delete...", el(session.text("Files---My-files---bdd-files-{time}---sub-b tree node inside browse tree"))));
      await session.step(52, "Then \"Are you sure?\" dialog should contain text \"bdd-files-{time}/sub-b\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("bdd-files-{time}/sub-b")));
      await session.step(53, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(54, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(55, "And Files---My-files---bdd-files-{time}---sub-b tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---sub-b tree node inside browse tree")), "absent"));
      await session.step(56, "And Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree")), "visible"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Create file... asks for a name and the tree shows the new file", async () => {
      await session.step(60, "When user picks \"Create file...\" from the context menu of Files---My-files---bdd-files-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Create file...", el(session.text("Files---My-files---bdd-files-{time} tree node inside browse tree"))));
      await session.step(61, "Then \"Create file\" dialog should be visible", () => shouldBe(page, el("\"Create file\" dialog"), "visible"));
      await session.step(62, "When user enters \"notes-a.txt\" into Name input in \"Create file\" dialog", () => enterInto(page, "notes-a.txt", el("Name input in \"Create file\" dialog")));
      await session.step(63, "And user clicks on OK button in \"Create file\" dialog", () => clickOn(page, el("OK button in \"Create file\" dialog")));
      await session.step(64, "Then the \"Create file\" dialog should close", () => dialogCloses(page, "Create file"));
      await session.step(65, "And Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree")), "visible"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Rename... shows the file's name and the tree shows the new one", async () => {
      await session.step(69, "When user picks \"Rename...\" from the context menu of Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree"))));
      await session.step(70, "Then \"File name\" input in Rename dialog should have value \"notes-a.txt\"", () => shouldHaveValue(page, el("\"File name\" input in Rename dialog"), "notes-a.txt"));
      await session.step(71, "When user enters \"notes-b.txt\" into \"File name\" input in Rename dialog", () => enterInto(page, "notes-b.txt", el("\"File name\" input in Rename dialog")));
      await session.step(72, "And user clicks on OK button in Rename dialog", () => clickOn(page, el("OK button in Rename dialog")));
      await session.step(73, "Then the \"Rename\" dialog should close", () => dialogCloses(page, "Rename"));
      await session.step(74, "And Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree")), "visible"));
      await session.step(75, "And Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree")), "absent"));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Delete... asks for confirmation and the file leaves the tree", async () => {
      await session.step(79, "When user picks \"Delete...\" from the context menu of Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree", () => pickFromContextMenu(page, "Delete...", el(session.text("Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree"))));
      await session.step(80, "Then \"Are you sure?\" dialog should contain text \"notes-b.txt\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "notes-b.txt"));
      await session.step(81, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(82, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(83, "And Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree")), "absent"));
      await session.step(84, "And Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree")), "visible"));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
