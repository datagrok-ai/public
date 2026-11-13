/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/save-to-space.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [file.menu.save.tables-as-project, views.space]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '../../bindings/identifiers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {projectInSpace, projectLinkedFromSpace, projectNotInSpace} from '../../bindings/spaces.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, pressKey, selectIn, shouldBe, shouldContainText, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, dialogCloses, noProjectOnServer, noSpaceOnServer, openDataset, openSaveDialog, saveFromDialog, simpleModeOff, spaceOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Space row of the Save project dialog puts a project into a space", () => {
  const session = feature(test, "features/spaces/save-to-space.feature", import.meta.url);
  test("The Space row of the Save project dialog puts a project into a space", {tag: ["@journey", "@realizes:file.menu.save.tables-as-project", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And simple mode is off", () => simpleModeOff(page));
    await session.step(9, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(10, "And a space named \"BDDSaveSpace{time}\" is on the server", () => spaceOnServer(page, session.text("BDDSaveSpace{time}")));
    await session.step(11, "And a space named \"BDDSaveLink{time}\" is on the server", () => spaceOnServer(page, session.text("BDDSaveLink{time}")));
    await session.step(12, "And a space named \"BDDSaveMove{time}\" is on the server", () => spaceOnServer(page, session.text("BDDSaveMove{time}")));
    await session.step(13, "And no project named \"BDDSpaceProj{time}, BDDNestedProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDSpaceProj{time}, BDDNestedProj{time}")));
    await session.step(14, "And no space named \"BDDSaveInner{time}, BDDSaveTop{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDSaveInner{time}, BDDSaveTop{time}")));
    await run.scenario("The Space row says what an empty choice means", async () => {
      await session.step(17, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(18, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
      await session.step(19, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(20, "And space input in \"Save project\" dialog should contain text \"Save to a space (optional)\"", () => shouldContainText(page, el("space input in \"Save project\" dialog"), "Save to a space (optional)"));
      await session.step(21, "When user hovers over space input in \"Save project\" dialog", () => hoverOver(page, el("space input in \"Save project\" dialog")));
      await session.step(22, "Then tooltip should contain text \"Empty: saved to your own projects\"", () => shouldContainText(page, el("tooltip"), "Empty: saved to your own projects"));
    });
    await run.scenario("A new project picked into a space is saved straight there", async () => {
      await session.step(25, "When user enters \"BDDSpaceProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDSpaceProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(26, "And user clicks on space input in \"Save project\" dialog", () => clickOn(page, el("space input in \"Save project\" dialog")));
      await session.step(27, "And user clicks on \"BDDSaveSpace{time}\" tree node", () => clickOn(page, el(session.text("\"BDDSaveSpace{time}\" tree node"))));
      await session.step(28, "Then space input in \"Save project\" dialog should contain text \"BDDSaveSpace{time}\"", () => shouldContainText(page, el("space input in \"Save project\" dialog"), session.text("BDDSaveSpace{time}")));
      await session.step(29, "And \"Save project\" dialog should contain text \"Shared with whoever can see \\\"BDDSaveSpace{time}\\\".\"", () => shouldContainText(page, el("\"Save project\" dialog"), session.text("Shared with whoever can see \"BDDSaveSpace{time}\".")));
      await session.step(30, "And \"Move\" choice chip in \"Save project\" dialog should be hidden", () => shouldBe(page, el("\"Move\" choice chip in \"Save project\" dialog"), "hidden"));
      await session.step(31, "When user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
      await session.step(32, "Then the project named \"BDDSpaceProj{time}\" should be in the space \"BDDSaveSpace{time}\" on the server", () => projectInSpace(page, session.text("BDDSpaceProj{time}"), session.text("BDDSaveSpace{time}")));
    });
    await run.scenario("Re-saved with Link, the project stays and the other space links it", async () => {
      await session.step(35, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
      await session.step(36, "Then space input in \"Save project\" dialog should contain text \"BDDSaveSpace{time}\"", () => shouldContainText(page, el("space input in \"Save project\" dialog"), session.text("BDDSaveSpace{time}")));
      await session.step(37, "And \"Save project\" dialog should contain text \"Lives in \\\"BDDSaveSpace{time}\\\".\"", () => shouldContainText(page, el("\"Save project\" dialog"), session.text("Lives in \"BDDSaveSpace{time}\".")));
      await session.step(38, "When user clicks on space input in \"Save project\" dialog", () => clickOn(page, el("space input in \"Save project\" dialog")));
      await session.step(39, "And user clicks on \"BDDSaveLink{time}\" tree node", () => clickOn(page, el(session.text("\"BDDSaveLink{time}\" tree node"))));
      await session.step(40, "Then \"Save project\" dialog should contain text \"Moves from \\\"BDDSaveSpace{time}\\\" to \\\"BDDSaveLink{time}\\\".\"", () => shouldContainText(page, el("\"Save project\" dialog"), session.text("Moves from \"BDDSaveSpace{time}\" to \"BDDSaveLink{time}\".")));
      await session.step(41, "When user clicks on \"Move\" choice chip in \"Save project\" dialog", () => clickOn(page, el("\"Move\" choice chip in \"Save project\" dialog")));
      await session.step(42, "And user clicks on \"Link\" menu item", () => clickOn(page, el("\"Link\" menu item")));
      await session.step(43, "Then \"Save project\" dialog should contain text \"Stays in \\\"BDDSaveSpace{time}\\\"; \\\"BDDSaveLink{time}\\\" gets a link to it.\"", () => shouldContainText(page, el("\"Save project\" dialog"), session.text("Stays in \"BDDSaveSpace{time}\"; \"BDDSaveLink{time}\" gets a link to it.")));
      await session.step(44, "When user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
      await session.step(45, "Then the project named \"BDDSpaceProj{time}\" should be in the space \"BDDSaveSpace{time}\" on the server", () => projectInSpace(page, session.text("BDDSpaceProj{time}"), session.text("BDDSaveSpace{time}")));
      await session.step(46, "And the project named \"BDDSpaceProj{time}\" should be linked from the space \"BDDSaveLink{time}\" on the server", () => projectLinkedFromSpace(page, session.text("BDDSpaceProj{time}"), session.text("BDDSaveLink{time}")));
    });
    await run.scenario("Moved into the space that links it, the project makes the link its home", async () => {
      await session.step(49, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
      await session.step(50, "And user clicks on space input in \"Save project\" dialog", () => clickOn(page, el("space input in \"Save project\" dialog")));
      await session.step(51, "And user clicks on \"BDDSaveLink{time}\" tree node", () => clickOn(page, el(session.text("\"BDDSaveLink{time}\" tree node"))));
      await session.step(52, "Then \"Save project\" dialog should contain text \"Moves from \\\"BDDSaveSpace{time}\\\" to \\\"BDDSaveLink{time}\\\".\"", () => shouldContainText(page, el("\"Save project\" dialog"), session.text("Moves from \"BDDSaveSpace{time}\" to \"BDDSaveLink{time}\".")));
      await session.step(53, "And \"Save project\" dialog should not contain text \"already has a project named\"", () => shouldNotContainText(page, el("\"Save project\" dialog"), "already has a project named"));
      await session.step(54, "And OK button in \"Save project\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Save project\" dialog"), "enabled"));
      await session.step(55, "When user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
      await session.step(56, "Then the project named \"BDDSpaceProj{time}\" should be in the space \"BDDSaveLink{time}\" on the server", () => projectInSpace(page, session.text("BDDSpaceProj{time}"), session.text("BDDSaveLink{time}")));
      await session.step(57, "And the project named \"BDDSpaceProj{time}\" should not be in the space \"BDDSaveSpace{time}\" on the server", () => projectNotInSpace(page, session.text("BDDSpaceProj{time}"), session.text("BDDSaveSpace{time}")));
    });
    await run.scenario("Re-saved with Move, the project goes to the other space", async () => {
      await session.step(60, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
      await session.step(61, "And user clicks on space input in \"Save project\" dialog", () => clickOn(page, el("space input in \"Save project\" dialog")));
      await session.step(62, "And user clicks on \"BDDSaveMove{time}\" tree node", () => clickOn(page, el(session.text("\"BDDSaveMove{time}\" tree node"))));
      await session.step(63, "Then \"Save project\" dialog should contain text \"Moves from \\\"BDDSaveLink{time}\\\" to \\\"BDDSaveMove{time}\\\".\"", () => shouldContainText(page, el("\"Save project\" dialog"), session.text("Moves from \"BDDSaveLink{time}\" to \"BDDSaveMove{time}\".")));
      await session.step(64, "When user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
      await session.step(65, "Then the project named \"BDDSpaceProj{time}\" should be in the space \"BDDSaveMove{time}\" on the server", () => projectInSpace(page, session.text("BDDSpaceProj{time}"), session.text("BDDSaveMove{time}")));
      await session.step(66, "And the project named \"BDDSpaceProj{time}\" should not be in the space \"BDDSaveLink{time}\" on the server", () => projectNotInSpace(page, session.text("BDDSpaceProj{time}"), session.text("BDDSaveLink{time}")));
    });
    await run.scenario("A space that holds a project of the name holds OK off", async () => {
      await session.step(69, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
      await session.step(70, "And user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
      await session.step(71, "And user enters \"BDDSpaceProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDSpaceProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(72, "And user clicks on space input in \"Save project\" dialog", () => clickOn(page, el("space input in \"Save project\" dialog")));
      await session.step(73, "And user clicks on \"BDDSaveMove{time}\" tree node", () => clickOn(page, el(session.text("\"BDDSaveMove{time}\" tree node"))));
      await session.step(74, "Then \"Save project\" dialog should contain text \"\\\"BDDSaveMove{time}\\\" already has a project named \\\"BDDSpaceProj{time}\\\"\"", () => shouldContainText(page, el("\"Save project\" dialog"), session.text("\"BDDSaveMove{time}\" already has a project named \"BDDSpaceProj{time}\"")));
      await session.step(75, "And OK button in \"Save project\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Save project\" dialog"), "disabled"));
      await session.step(76, "When user enters \"BDDSpaceProjCopy{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDSpaceProjCopy{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(77, "Then \"Save project\" dialog should not contain text \"already has a project named\"", () => shouldNotContainText(page, el("\"Save project\" dialog"), "already has a project named"));
      await session.step(78, "And OK button in \"Save project\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Save project\" dialog"), "enabled"));
      await session.step(79, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(80, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    });
    await run.scenario("Spaces made from the \"+\" buttons nest, and the project is saved into the inner one", async () => {
      await session.step(83, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
      await session.step(84, "And user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
      await session.step(85, "And user enters \"BDDNestedProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDNestedProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(86, "And user clicks on space input in \"Save project\" dialog", () => clickOn(page, el("space input in \"Save project\" dialog")));
      await session.step(87, "And user clicks on \"New top-level space\" button", () => clickOn(page, el("\"New top-level space\" button")));
      await session.step(88, "And user types \"BDDSaveTop{time}\" into \"Rename New space\" text input", () => typeInto(page, session.text("BDDSaveTop{time}"), el("\"Rename New space\" text input")));
      await session.step(89, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(90, "Then space input in \"Save project\" dialog should contain text \"BDDSaveTop{time}\"", () => shouldContainText(page, el("space input in \"Save project\" dialog"), session.text("BDDSaveTop{time}")));
      await session.step(91, "When user clicks on \"New space inside BDDSaveTop{time}\" button", () => clickOn(page, el(session.text("\"New space inside BDDSaveTop{time}\" button"))));
      await session.step(92, "And user types \"BDDSaveInner{time}\" into \"Rename New space\" text input", () => typeInto(page, session.text("BDDSaveInner{time}"), el("\"Rename New space\" text input")));
      await session.step(93, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(94, "Then space input in \"Save project\" dialog should contain text \"BDDSaveTop{time} / BDDSaveInner{time}\"", () => shouldContainText(page, el("space input in \"Save project\" dialog"), session.text("BDDSaveTop{time} / BDDSaveInner{time}")));
      await session.step(95, "When user clicks on space input in \"Save project\" dialog", () => clickOn(page, el("space input in \"Save project\" dialog")));
      await session.step(96, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
      await session.step(97, "Then the project named \"BDDNestedProj{time}\" should be in the space \"BDDSaveInner{time}\" on the server", () => projectInSpace(page, session.text("BDDNestedProj{time}"), session.text("BDDSaveInner{time}")));
    });
    run.finish();
  });
});
