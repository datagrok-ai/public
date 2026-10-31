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

test.describe("Files: a folder and a file created, and a file renamed and deleted, from the Browse tree", () => {
  const session = feature(test, "features/browse/files-lifecycle.feature", import.meta.url);
  test("Files: a folder and a file created, and a file renamed and deleted, from the Browse tree", {tag: ["@journey", "@browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the \"fixtures/cars-small.csv\" file of the project is in the user's files as \"bdd-files-{time}/seed.csv\"", () => fixtureInUsersFiles(page, "fixtures/cars-small.csv", session.text("bdd-files-{time}/seed.csv")));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(26, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(27, "And Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
    await session.step(28, "And Files---My-files---bdd-files-{time} tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("Files---My-files---bdd-files-{time} tree node inside browse tree"))));
    await run.scenario("Create folder... and Create file... ask for a name, and the tree shows the new folder and file", async () => {
      await session.step(31, "Then Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree")), "visible"));
      await session.step(32, "When user picks \"Create folder...\" from the context menu of Files---My-files---bdd-files-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Create folder...", el(session.text("Files---My-files---bdd-files-{time} tree node inside browse tree"))));
      await session.step(33, "Then \"Create folder\" dialog should be visible", () => shouldBe(page, el("\"Create folder\" dialog"), "visible"));
      await session.step(34, "When user enters \"sub-a\" into Name input in \"Create folder\" dialog", () => enterInto(page, "sub-a", el("Name input in \"Create folder\" dialog")));
      await session.step(35, "And user clicks on OK button in \"Create folder\" dialog", () => clickOn(page, el("OK button in \"Create folder\" dialog")));
      await session.step(36, "Then the \"Create folder\" dialog should close", () => dialogCloses(page, "Create folder"));
      await session.step(37, "And Files---My-files---bdd-files-{time}---sub-a tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---sub-a tree node inside browse tree")), "visible"));
      await session.step(38, "When user picks \"Create file...\" from the context menu of Files---My-files---bdd-files-{time} tree node inside browse tree", () => pickFromContextMenu(page, "Create file...", el(session.text("Files---My-files---bdd-files-{time} tree node inside browse tree"))));
      await session.step(39, "Then \"Create file\" dialog should be visible", () => shouldBe(page, el("\"Create file\" dialog"), "visible"));
      await session.step(40, "When user enters \"notes-a.txt\" into Name input in \"Create file\" dialog", () => enterInto(page, "notes-a.txt", el("Name input in \"Create file\" dialog")));
      await session.step(41, "And user clicks on OK button in \"Create file\" dialog", () => clickOn(page, el("OK button in \"Create file\" dialog")));
      await session.step(42, "Then the \"Create file\" dialog should close", () => dialogCloses(page, "Create file"));
      await session.step(43, "And Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree")), "visible"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Rename... shows the file's name and renames it, and Delete... asks for confirmation and removes it", async () => {
      await session.step(47, "When user picks \"Rename...\" from the context menu of Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree"))));
      await session.step(48, "Then \"File name\" input in Rename dialog should have value \"notes-a.txt\"", () => shouldHaveValue(page, el("\"File name\" input in Rename dialog"), "notes-a.txt"));
      await session.step(49, "When user enters \"notes-b.txt\" into \"File name\" input in Rename dialog", () => enterInto(page, "notes-b.txt", el("\"File name\" input in Rename dialog")));
      await session.step(50, "And user clicks on OK button in Rename dialog", () => clickOn(page, el("OK button in Rename dialog")));
      await session.step(51, "Then the \"Rename\" dialog should close", () => dialogCloses(page, "Rename"));
      await session.step(52, "And Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree")), "visible"));
      await session.step(53, "And Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---notes-a.txt tree node inside browse tree")), "absent"));
      await session.step(54, "When user picks \"Delete...\" from the context menu of Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree", () => pickFromContextMenu(page, "Delete...", el(session.text("Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree"))));
      await session.step(55, "Then \"Are you sure?\" dialog should contain text \"notes-b.txt\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "notes-b.txt"));
      await session.step(56, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(57, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(58, "And Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---notes-b.txt tree node inside browse tree")), "absent"));
      await session.step(59, "And Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("Files---My-files---bdd-files-{time}---seed.csv tree node inside browse tree")), "visible"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
