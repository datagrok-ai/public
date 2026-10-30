/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/campaign/progress-tracker.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addStage, closeViewByName, hitDesignCampaign, noLockHeld, openByLink, savedCell, savedVersion} from '../../bindings/campaign.js';
import {inside} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, shouldBe, shouldOffer} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDatasetRows, switchView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea, infoBalloonText, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragAreaOntoWidget} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The progress tracker of a Hit Design campaign", () => {
  const session = feature(test, "features/campaign/progress-tracker.feature", import.meta.url);
  test("The tracker covers the campaign view with a lane per stage, in the template's order, empty ones included", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the Hit Design campaign \"BDD Tracker\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Tracker", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]);
    await session.step(19, "When user opens the campaign \"BDD Tracker\" by its link", () => openByLink(page, "BDD Tracker"));
    await session.step(20, "And user clicks on \"Progress tracker\" button in toolbar", () => clickOn(page, el("\"Progress tracker\" button in toolbar")));
    await session.step(21, "Then \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "visible"));
    await session.step(24, "Then grid should be inside \"Progress Tracker\" dialog", () => inside(page, el("grid"), el("\"Progress Tracker\" dialog")));
    await session.step(25, "And the \"lane names\" reading of \"Tile Viewer\" viewer should be \"Idea, Design, Make\"", () => readingReads(page, "lane names", el("\"Tile Viewer\" viewer"), "Idea, Design, Make"));
    await session.step(26, "And the 'tiles in lane \"Make\"' reading of \"Tile Viewer\" viewer should be 0", () => readingIs(page, "tiles in lane \"Make\"", el("\"Tile Viewer\" viewer"), 0));
    await session.step(27, "And the \"lane of row 1\" reading of \"Tile Viewer\" viewer should be \"Idea\"", () => readingReads(page, "lane of row 1", el("\"Tile Viewer\" viewer"), "Idea"));
    await session.step(28, "And the \"lane of row 2\" reading of \"Tile Viewer\" viewer should be \"Design\"", () => readingReads(page, "lane of row 2", el("\"Tile Viewer\" viewer"), "Design"));
    await session.step(29, "And CLOSE button in \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("CLOSE button in \"Progress Tracker\" dialog"), "visible"));
    await session.step(30, "And \"Modify Stages\" button in \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("\"Modify Stages\" button in \"Progress Tracker\" dialog"), "visible"));
    await session.step(31, "And CANCEL button in \"Progress Tracker\" dialog should be absent", () => shouldBe(page, el("CANCEL button in \"Progress Tracker\" dialog"), "absent"));
    await session.step(32, "And OK button in \"Progress Tracker\" dialog should be absent", () => shouldBe(page, el("OK button in \"Progress Tracker\" dialog"), "absent"));
  });
  test("A card dragged to another lane moves its row there, and CLOSE saves it with a balloon", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the Hit Design campaign \"BDD Tracker\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Tracker", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]);
    await session.step(19, "When user opens the campaign \"BDD Tracker\" by its link", () => openByLink(page, "BDD Tracker"));
    await session.step(20, "And user clicks on \"Progress tracker\" button in toolbar", () => clickOn(page, el("\"Progress tracker\" button in toolbar")));
    await session.step(21, "Then \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "visible"));
    await session.step(35, "When user drags the \"tile of row 1\" area of \"Tile Viewer\" viewer onto the 'lane content \"Make\"' area of \"Tile Viewer\" viewer", () => dragAreaOntoWidget(page, "tile of row 1", el("\"Tile Viewer\" viewer"), "lane content \"Make\"", el("\"Tile Viewer\" viewer")));
    await session.step(36, "Then the \"lane of row 1\" reading of \"Tile Viewer\" viewer should be \"Make\"", () => readingReads(page, "lane of row 1", el("\"Tile Viewer\" viewer"), "Make"));
    await session.step(37, "And the value of \"Stage\" column in row 1 should be \"Make\"", () => valueInRow(page, "Stage", 1, "Make"));
    await session.step(38, "When user clicks on CLOSE button in \"Progress Tracker\" dialog", () => clickOn(page, el("CLOSE button in \"Progress Tracker\" dialog")));
    await session.step(39, "Then \"Progress Tracker\" dialog should be absent", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "absent"));
    await session.step(40, "And an info balloon containing \"Campaign saved successfully\" should have been shown", () => infoBalloonText(page, "Campaign saved successfully"));
    await session.step(41, "And the campaign \"BDD Tracker\" should be saved with \"Make\" in row 1 of the \"Stage\" column", () => savedCell(page, "BDD Tracker", "Make", 1, "Stage"));
  });
  test("Modify Stages adds a stage, saves the campaign and reopens the tracker with the new lane", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the Hit Design campaign \"BDD Tracker\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Tracker", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]);
    await session.step(19, "When user opens the campaign \"BDD Tracker\" by its link", () => openByLink(page, "BDD Tracker"));
    await session.step(20, "And user clicks on \"Progress tracker\" button in toolbar", () => clickOn(page, el("\"Progress tracker\" button in toolbar")));
    await session.step(21, "Then \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "visible"));
    await session.step(45, "When user clicks on \"Modify Stages\" button in \"Progress Tracker\" dialog", () => clickOn(page, el("\"Modify Stages\" button in \"Progress Tracker\" dialog")));
    await session.step(46, "Then \"Modify Stages\" dialog should be visible", () => shouldBe(page, el("\"Modify Stages\" dialog"), "visible"));
    await session.step(47, "And CANCEL button in \"Modify Stages\" dialog should be enabled", () => shouldBe(page, el("CANCEL button in \"Modify Stages\" dialog"), "enabled"));
    await session.step(48, "When user adds the stage \"Test\" in \"Modify Stages\" dialog", () => addStage(page, "Test", el("\"Modify Stages\" dialog")));
    await session.step(49, "And user clicks on OK button in \"Modify Stages\" dialog", () => clickOn(page, el("OK button in \"Modify Stages\" dialog")));
    await session.step(50, "Then \"Modify Stages\" dialog should be absent", () => shouldBe(page, el("\"Modify Stages\" dialog"), "absent"));
    await session.step(51, "And an info balloon containing \"Campaign saved successfully\" should have been shown", () => infoBalloonText(page, "Campaign saved successfully"));
    await session.step(52, "And \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "visible"));
    await session.step(53, "And the \"lane names\" reading of \"Tile Viewer\" viewer should be \"Idea, Design, Make, Test\"", () => readingReads(page, "lane names", el("\"Tile Viewer\" viewer"), "Idea, Design, Make, Test"));
    await session.step(54, "When user clicks on CLOSE button in \"Progress Tracker\" dialog", () => clickOn(page, el("CLOSE button in \"Progress Tracker\" dialog")));
    await session.step(55, "Then \"Progress Tracker\" dialog should be absent", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "absent"));
    await session.step(56, "When user double-clicks on the \"cell 1 of Stage\" area of grid", () => doubleClickArea(page, "cell 1 of Stage", el("grid")));
    await session.step(57, "Then cell editor should offer \"Idea, Design, Make, Test\"", () => shouldOffer(page, el("cell editor"), "Idea, Design, Make, Test"));
  });
  test("The tracker hides while another view is current and comes back with the campaign view", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the Hit Design campaign \"BDD Tracker\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Tracker", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]);
    await session.step(19, "When user opens the campaign \"BDD Tracker\" by its link", () => openByLink(page, "BDD Tracker"));
    await session.step(20, "And user clicks on \"Progress tracker\" button in toolbar", () => clickOn(page, el("\"Progress tracker\" button in toolbar")));
    await session.step(21, "Then \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "visible"));
    await session.step(60, "Given user opens demog dataset keeping the first 10 rows", () => openDatasetRows(page, ds("demog"), 10));
    await session.step(61, "Then \"Progress Tracker\" dialog should be hidden", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "hidden"));
    await session.step(62, "When user switches to the \"BDD Tracker\" view", () => switchView(page, "BDD Tracker"));
    await session.step(63, "Then \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "visible"));
    await session.step(64, "And the \"lane names\" reading of \"Tile Viewer\" viewer should be \"Idea, Design, Make\"", () => readingReads(page, "lane names", el("\"Tile Viewer\" viewer"), "Idea, Design, Make"));
  });
  test("Closing the campaign view while its tracker is hidden takes the tracker away", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the Hit Design campaign \"BDD Tracker\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Tracker", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]);
    await session.step(19, "When user opens the campaign \"BDD Tracker\" by its link", () => openByLink(page, "BDD Tracker"));
    await session.step(20, "And user clicks on \"Progress tracker\" button in toolbar", () => clickOn(page, el("\"Progress tracker\" button in toolbar")));
    await session.step(21, "Then \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "visible"));
    await session.step(67, "Given user opens demog dataset keeping the first 10 rows", () => openDatasetRows(page, ds("demog"), 10));
    await session.step(68, "Then \"Progress Tracker\" dialog should be hidden", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "hidden"));
    await session.step(69, "When user closes the \"BDD Tracker\" view", () => closeViewByName(page, "BDD Tracker"));
    await session.step(70, "Then \"Progress Tracker\" dialog should be absent", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "absent"));
  });
  test("Escape and the close button keep the moves and save them, as CLOSE does", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the Hit Design campaign \"BDD Tracker\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Tracker", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]);
    await session.step(19, "When user opens the campaign \"BDD Tracker\" by its link", () => openByLink(page, "BDD Tracker"));
    await session.step(20, "And user clicks on \"Progress tracker\" button in toolbar", () => clickOn(page, el("\"Progress tracker\" button in toolbar")));
    await session.step(21, "Then \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "visible"));
    await session.step(73, "When user drags the \"tile of row 1\" area of \"Tile Viewer\" viewer onto the 'lane content \"Make\"' area of \"Tile Viewer\" viewer", () => dragAreaOntoWidget(page, "tile of row 1", el("\"Tile Viewer\" viewer"), "lane content \"Make\"", el("\"Tile Viewer\" viewer")));
    await session.step(74, "Then the \"lane of row 1\" reading of \"Tile Viewer\" viewer should be \"Make\"", () => readingReads(page, "lane of row 1", el("\"Tile Viewer\" viewer"), "Make"));
    await session.step(75, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(76, "Then \"Progress Tracker\" dialog should be absent", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "absent"));
    await session.step(77, "And the value of \"Stage\" column in row 1 should be \"Make\"", () => valueInRow(page, "Stage", 1, "Make"));
    await session.step(78, "And the campaign \"BDD Tracker\" should be at version 2 on the server", () => savedVersion(page, "BDD Tracker", 2));
    await session.step(79, "And the campaign \"BDD Tracker\" should be saved with \"Make\" in row 1 of the \"Stage\" column", () => savedCell(page, "BDD Tracker", "Make", 1, "Stage"));
    await session.step(80, "When user clicks on \"Progress tracker\" button in toolbar", () => clickOn(page, el("\"Progress tracker\" button in toolbar")));
    await session.step(81, "Then \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "visible"));
    await session.step(82, "When user drags the \"tile of row 2\" area of \"Tile Viewer\" viewer onto the 'lane content \"Idea\"' area of \"Tile Viewer\" viewer", () => dragAreaOntoWidget(page, "tile of row 2", el("\"Tile Viewer\" viewer"), "lane content \"Idea\"", el("\"Tile Viewer\" viewer")));
    await session.step(83, "Then the \"lane of row 2\" reading of \"Tile Viewer\" viewer should be \"Idea\"", () => readingReads(page, "lane of row 2", el("\"Tile Viewer\" viewer"), "Idea"));
    await session.step(84, "When user clicks on close button of \"Progress Tracker\" dialog", () => clickOn(page, el("close button of \"Progress Tracker\" dialog")));
    await session.step(85, "Then \"Progress Tracker\" dialog should be absent", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "absent"));
    await session.step(86, "And the value of \"Stage\" column in row 2 should be \"Idea\"", () => valueInRow(page, "Stage", 2, "Idea"));
    await session.step(87, "And the campaign \"BDD Tracker\" should be at version 3 on the server", () => savedVersion(page, "BDD Tracker", 3));
    await session.step(88, "And the campaign \"BDD Tracker\" should be saved with \"Idea\" in row 2 of the \"Stage\" column", () => savedCell(page, "BDD Tracker", "Idea", 2, "Stage"));
    await session.step(89, "And the campaign \"BDD Tracker\" should hold no lock on the server", () => noLockHeld(page, "BDD Tracker"));
  });
});
