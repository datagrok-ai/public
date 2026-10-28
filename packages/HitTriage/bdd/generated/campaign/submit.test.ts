/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/campaign/submit.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {hitDesignCampaignWithMolecule, openByLink, savedName, savedStatus, structureInRow, suggests, summaryIs, takeSuggestion} from '../../bindings/campaign.js';
import {clearOf} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, downloadThrough, downloadedContains, fileDownloaded, pressKeyIn, selectIn, shouldBe, shouldContainText, shouldHaveValue, typeInto, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {shouldHaveError, shouldHaveNoError} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {errorBalloonText, infoBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Submit dialog of a Hit Design campaign", () => {
  const session = feature(test, "features/campaign/submit.feature", import.meta.url);
  test("The summary lists the campaign's keys in order, and CANCEL closes the dialog", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the Hit Design campaign \"BDD Submit\" with the stages \"Idea, Design\", the molecule field \"Reference compound\" set to \"c1ccc2[nH]ccc2c1\" and a table of:", () => hitDesignCampaignWithMolecule(page, "BDD Submit", "Idea, Design", "Reference compound", "c1ccc2[nH]ccc2c1", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]);
    await session.step(25, "When user opens the campaign \"BDD Submit\" by its link", () => openByLink(page, "BDD Submit"));
    await session.step(26, "And user clicks on \"Submit...\" button in toolbar", () => clickOn(page, el("\"Submit...\" button in toolbar")));
    await session.step(27, "Then \"Submit\" dialog should be visible", () => shouldBe(page, el("\"Submit\" dialog"), "visible"));
    await session.step(30, "Then the summary in \"Submit\" dialog should be:", () => summaryIs(page, el("\"Submit\" dialog"), [["Name","BDD Submit"],["Code","BDDBDDSUBMIT-1"],["Template","BDD Submit"],["File Path","System.AppData/HitTriage/Hit Design/campaigns…"],["Create Date","2026/10/15"],["Chemist","BDD chemist"],["Reference compound","*"],["Number of Molecules","2"],["Enrichment Methods","*"]]), [["Name","BDD Submit"],["Code","BDDBDDSUBMIT-1"],["Template","BDD Submit"],["File Path","System.AppData/HitTriage/Hit Design/campaigns…"],["Create Date","2026/10/15"],["Chemist","BDD chemist"],["Reference compound","*"],["Number of Molecules","2"],["Enrichment Methods","*"]]);
    await session.step(40, "And \"Status\" input in \"Submit\" dialog should have value \"In Progress\"", () => shouldHaveValue(page, el("\"Status\" input in \"Submit\" dialog"), "In Progress"));
    await session.step(41, "And \"Submit function\" input in \"Submit\" dialog should have value \"Demo File Submit\"", () => shouldHaveValue(page, el("\"Submit function\" input in \"Submit\" dialog"), "Demo File Submit"));
    await session.step(42, "And CANCEL button in \"Submit\" dialog should be enabled", () => shouldBe(page, el("CANCEL button in \"Submit\" dialog"), "enabled"));
    await session.step(43, "When user clicks on CANCEL button in \"Submit\" dialog", () => clickOn(page, el("CANCEL button in \"Submit\" dialog")));
    await session.step(44, "Then \"Submit\" dialog should be absent", () => shouldBe(page, el("\"Submit\" dialog"), "absent"));
  });
  test("The status suggests the statuses campaigns already have, and Save stores it with a balloon", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the Hit Design campaign \"BDD Submit\" with the stages \"Idea, Design\", the molecule field \"Reference compound\" set to \"c1ccc2[nH]ccc2c1\" and a table of:", () => hitDesignCampaignWithMolecule(page, "BDD Submit", "Idea, Design", "Reference compound", "c1ccc2[nH]ccc2c1", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]);
    await session.step(25, "When user opens the campaign \"BDD Submit\" by its link", () => openByLink(page, "BDD Submit"));
    await session.step(26, "And user clicks on \"Submit...\" button in toolbar", () => clickOn(page, el("\"Submit...\" button in toolbar")));
    await session.step(27, "Then \"Submit\" dialog should be visible", () => shouldBe(page, el("\"Submit\" dialog"), "visible"));
    await session.step(47, "When user clears \"Status\" input in \"Submit\" dialog", () => clearField(page, el("\"Status\" input in \"Submit\" dialog")));
    await session.step(48, "And user types \"Progr\" into \"Status\" input in \"Submit\" dialog", () => typeInto(page, "Progr", el("\"Status\" input in \"Submit\" dialog")));
    await session.step(49, "Then \"Status\" input in \"Submit\" dialog should suggest \"In Progress\"", () => suggests(page, el("\"Status\" input in \"Submit\" dialog"), "In Progress"));
    await session.step(50, "When user takes the suggestion \"In Progress\" of \"Status\" input in \"Submit\" dialog", () => takeSuggestion(page, "In Progress", el("\"Status\" input in \"Submit\" dialog")));
    await session.step(51, "Then \"Status\" input in \"Submit\" dialog should have value \"In Progress\"", () => shouldHaveValue(page, el("\"Status\" input in \"Submit\" dialog"), "In Progress"));
    await session.step(52, "When user clears \"Status\" input in \"Submit\" dialog", () => clearField(page, el("\"Status\" input in \"Submit\" dialog")));
    await session.step(53, "And user types \"BDD reviewed\" into \"Status\" input in \"Submit\" dialog", () => typeInto(page, "BDD reviewed", el("\"Status\" input in \"Submit\" dialog")));
    await session.step(54, "And user clicks on Save button in \"Submit\" dialog", () => clickOn(page, el("Save button in \"Submit\" dialog")));
    await session.step(55, "Then \"Submit\" dialog should be absent", () => shouldBe(page, el("\"Submit\" dialog"), "absent"));
    await session.step(56, "And an info balloon containing \"Campaign saved successfully\" should have been shown", () => infoBalloonText(page, "Campaign saved successfully"));
    await session.step(57, "And the campaign \"BDD Submit\" should be saved with the status \"BDD reviewed\"", () => savedStatus(page, "BDD Submit", "BDD reviewed"));
  });
  test("Submit calls the demo submit function with the campaign's rows, saves the status and says so", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the Hit Design campaign \"BDD Submit\" with the stages \"Idea, Design\", the molecule field \"Reference compound\" set to \"c1ccc2[nH]ccc2c1\" and a table of:", () => hitDesignCampaignWithMolecule(page, "BDD Submit", "Idea, Design", "Reference compound", "c1ccc2[nH]ccc2c1", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]);
    await session.step(25, "When user opens the campaign \"BDD Submit\" by its link", () => openByLink(page, "BDD Submit"));
    await session.step(26, "And user clicks on \"Submit...\" button in toolbar", () => clickOn(page, el("\"Submit...\" button in toolbar")));
    await session.step(27, "Then \"Submit\" dialog should be visible", () => shouldBe(page, el("\"Submit\" dialog"), "visible"));
    await session.step(60, "When user clears \"Status\" input in \"Submit\" dialog", () => clearField(page, el("\"Status\" input in \"Submit\" dialog")));
    await session.step(61, "And user types \"Submitted\" into \"Status\" input in \"Submit\" dialog", () => typeInto(page, "Submitted", el("\"Status\" input in \"Submit\" dialog")));
    await session.step(62, "And user clicks on Submit button in \"Submit\" dialog", () => clickOn(page, el("Submit button in \"Submit\" dialog")));
    await session.step(63, "Then \"Submit\" dialog should be absent", () => shouldBe(page, el("\"Submit\" dialog"), "absent"));
    await session.step(64, "And an info balloon containing \"Submitted successfully\" should have been shown", () => infoBalloonText(page, "Submitted successfully"));
    await session.step(65, "And an info balloon containing \"Campaign saved successfully\" should have been shown", () => infoBalloonText(page, "Campaign saved successfully"));
    await session.step(66, "And an info balloon containing \"Molecule\" should have been shown", () => infoBalloonText(page, "Molecule"));
    await session.step(67, "And the campaign \"BDD Submit\" should be saved with the status \"Submitted\"", () => savedStatus(page, "BDD Submit", "Submitted"));
  });
  test("Without a submit function Submit stays, disabled, saying why", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the Hit Design campaign \"BDD Submit\" with the stages \"Idea, Design\", the molecule field \"Reference compound\" set to \"c1ccc2[nH]ccc2c1\" and a table of:", () => hitDesignCampaignWithMolecule(page, "BDD Submit", "Idea, Design", "Reference compound", "c1ccc2[nH]ccc2c1", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]);
    await session.step(25, "When user opens the campaign \"BDD Submit\" by its link", () => openByLink(page, "BDD Submit"));
    await session.step(26, "And user clicks on \"Submit...\" button in toolbar", () => clickOn(page, el("\"Submit...\" button in toolbar")));
    await session.step(27, "Then \"Submit\" dialog should be visible", () => shouldBe(page, el("\"Submit\" dialog"), "visible"));
    await session.step(70, "When user selects \"\" in \"Submit function\" input in \"Submit\" dialog", () => selectIn(page, "", el("\"Submit function\" input in \"Submit\" dialog")));
    await session.step(71, "Then Save button in \"Submit\" dialog should be enabled", () => shouldBe(page, el("Save button in \"Submit\" dialog"), "enabled"));
    await session.step(72, "And Submit button in \"Submit\" dialog should be disabled", () => shouldBe(page, el("Submit button in \"Submit\" dialog"), "disabled"));
    await session.step(73, "And \"Submit\" dialog should contain text \"Choose a submit function\"", () => shouldContainText(page, el("\"Submit\" dialog"), "Choose a submit function"));
    await session.step(74, "When user selects \"Demo File Submit\" in \"Submit function\" input in \"Submit\" dialog", () => selectIn(page, "Demo File Submit", el("\"Submit function\" input in \"Submit\" dialog")));
    await session.step(75, "Then Submit button in \"Submit\" dialog should be enabled", () => shouldBe(page, el("Submit button in \"Submit\" dialog"), "enabled"));
  });
  test("Enrichment Methods names the template's compute by its friendly name", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the Hit Design campaign \"BDD Submit\" with the stages \"Idea, Design\", the molecule field \"Reference compound\" set to \"c1ccc2[nH]ccc2c1\" and a table of:", () => hitDesignCampaignWithMolecule(page, "BDD Submit", "Idea, Design", "Reference compound", "c1ccc2[nH]ccc2c1", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]);
    await session.step(25, "When user opens the campaign \"BDD Submit\" by its link", () => openByLink(page, "BDD Submit"));
    await session.step(26, "And user clicks on \"Submit...\" button in toolbar", () => clickOn(page, el("\"Submit...\" button in toolbar")));
    await session.step(27, "Then \"Submit\" dialog should be visible", () => shouldBe(page, el("\"Submit\" dialog"), "visible"));
    await session.step(78, "Then the summary in \"Submit\" dialog should be:", () => summaryIs(page, el("\"Submit\" dialog"), [["Name","*"],["Code","*"],["Template","*"],["File Path","*"],["Create Date","*"],["Chemist","*"],["Reference compound","*"],["Number of Molecules","*"],["Enrichment Methods","Chemical Properties"]]), [["Name","*"],["Code","*"],["Template","*"],["File Path","*"],["Create Date","*"],["Chemist","*"],["Reference compound","*"],["Number of Molecules","*"],["Enrichment Methods","Chemical Properties"]]);
  });
  test("A Molecule field shows its structure in the value column, on a row as tall as the picture", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the Hit Design campaign \"BDD Submit\" with the stages \"Idea, Design\", the molecule field \"Reference compound\" set to \"c1ccc2[nH]ccc2c1\" and a table of:", () => hitDesignCampaignWithMolecule(page, "BDD Submit", "Idea, Design", "Reference compound", "c1ccc2[nH]ccc2c1", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]);
    await session.step(25, "When user opens the campaign \"BDD Submit\" by its link", () => openByLink(page, "BDD Submit"));
    await session.step(26, "And user clicks on \"Submit...\" button in toolbar", () => clickOn(page, el("\"Submit...\" button in toolbar")));
    await session.step(27, "Then \"Submit\" dialog should be visible", () => shouldBe(page, el("\"Submit\" dialog"), "visible"));
    await session.step(90, "Then the structure of the \"Reference compound\" row in \"Submit\" dialog should sit in the value column of its row", () => structureInRow(page, "Reference compound", el("\"Submit\" dialog")));
    await session.step(91, "And \"Reference compound\" input in \"Submit\" dialog should not overlap \"Number of Molecules\" input in \"Submit\" dialog", () => clearOf(page, el("\"Reference compound\" input in \"Submit\" dialog"), el("\"Number of Molecules\" input in \"Submit\" dialog")));
  });
  test("The pencil opens the name; Escape and Cancel put it back, an empty name is refused under it", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the Hit Design campaign \"BDD Submit\" with the stages \"Idea, Design\", the molecule field \"Reference compound\" set to \"c1ccc2[nH]ccc2c1\" and a table of:", () => hitDesignCampaignWithMolecule(page, "BDD Submit", "Idea, Design", "Reference compound", "c1ccc2[nH]ccc2c1", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]);
    await session.step(25, "When user opens the campaign \"BDD Submit\" by its link", () => openByLink(page, "BDD Submit"));
    await session.step(26, "And user clicks on \"Submit...\" button in toolbar", () => clickOn(page, el("\"Submit...\" button in toolbar")));
    await session.step(27, "Then \"Submit\" dialog should be visible", () => shouldBe(page, el("\"Submit\" dialog"), "visible"));
    await session.step(95, "Then \"Name\" input in \"Submit\" dialog should be read-only", () => shouldBe(page, el("\"Name\" input in \"Submit\" dialog"), "read-only"));
    await session.step(96, "When user clicks on \"Edit Campaign Name\" icon in \"Submit\" dialog", () => clickOn(page, el("\"Edit Campaign Name\" icon in \"Submit\" dialog")));
    await session.step(97, "Then \"Name\" input in \"Submit\" dialog should be editable", () => shouldBe(page, el("\"Name\" input in \"Submit\" dialog"), "editable"));
    await session.step(98, "And \"Name\" input in \"Submit\" dialog should be focused", () => shouldBe(page, el("\"Name\" input in \"Submit\" dialog"), "focused"));
    await session.step(99, "When user types \"BDD Submit renamed\" into \"Name\" input in \"Submit\" dialog", () => typeInto(page, "BDD Submit renamed", el("\"Name\" input in \"Submit\" dialog")));
    await session.step(100, "And user presses Escape in \"Name\" input in \"Submit\" dialog", () => pressKeyIn(page, "Escape", el("\"Name\" input in \"Submit\" dialog")));
    await session.step(101, "Then \"Submit\" dialog should be visible", () => shouldBe(page, el("\"Submit\" dialog"), "visible"));
    await session.step(102, "And \"Name\" input in \"Submit\" dialog should have value \"BDD Submit\"", () => shouldHaveValue(page, el("\"Name\" input in \"Submit\" dialog"), "BDD Submit"));
    await session.step(103, "And \"Name\" input in \"Submit\" dialog should be read-only", () => shouldBe(page, el("\"Name\" input in \"Submit\" dialog"), "read-only"));
    await session.step(104, "When user clicks on \"Edit Campaign Name\" icon in \"Submit\" dialog", () => clickOn(page, el("\"Edit Campaign Name\" icon in \"Submit\" dialog")));
    await session.step(105, "And user clears \"Name\" input in \"Submit\" dialog", () => clearField(page, el("\"Name\" input in \"Submit\" dialog")));
    await session.step(106, "And user clicks on \"Save\" icon in \"Submit\" dialog", () => clickOn(page, el("\"Save\" icon in \"Submit\" dialog")));
    await session.step(107, "Then \"Name\" input in \"Submit\" dialog should have the error \"Campaign name can not be empty\"", () => shouldHaveError(page, el("\"Name\" input in \"Submit\" dialog"), "Campaign name can not be empty"));
    await session.step(108, "When user clicks on \"Cancel\" icon in \"Submit\" dialog", () => clickOn(page, el("\"Cancel\" icon in \"Submit\" dialog")));
    await session.step(109, "Then \"Name\" input in \"Submit\" dialog should have value \"BDD Submit\"", () => shouldHaveValue(page, el("\"Name\" input in \"Submit\" dialog"), "BDD Submit"));
    await session.step(110, "And \"Name\" input in \"Submit\" dialog should have no error", () => shouldHaveNoError(page, el("\"Name\" input in \"Submit\" dialog")));
  });
  test("Enter saves the new name with the campaign and leaves the dialog open", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the Hit Design campaign \"BDD Submit\" with the stages \"Idea, Design\", the molecule field \"Reference compound\" set to \"c1ccc2[nH]ccc2c1\" and a table of:", () => hitDesignCampaignWithMolecule(page, "BDD Submit", "Idea, Design", "Reference compound", "c1ccc2[nH]ccc2c1", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]);
    await session.step(25, "When user opens the campaign \"BDD Submit\" by its link", () => openByLink(page, "BDD Submit"));
    await session.step(26, "And user clicks on \"Submit...\" button in toolbar", () => clickOn(page, el("\"Submit...\" button in toolbar")));
    await session.step(27, "Then \"Submit\" dialog should be visible", () => shouldBe(page, el("\"Submit\" dialog"), "visible"));
    await session.step(113, "When user clicks on \"Edit Campaign Name\" icon in \"Submit\" dialog", () => clickOn(page, el("\"Edit Campaign Name\" icon in \"Submit\" dialog")));
    await session.step(114, "And user types \"BDD Submit renamed\" into \"Name\" input in \"Submit\" dialog", () => typeInto(page, "BDD Submit renamed", el("\"Name\" input in \"Submit\" dialog")));
    await session.step(115, "And user presses Enter in \"Name\" input in \"Submit\" dialog", () => pressKeyIn(page, "Enter", el("\"Name\" input in \"Submit\" dialog")));
    await session.step(116, "Then an info balloon containing \"Campaign saved successfully\" should have been shown", () => infoBalloonText(page, "Campaign saved successfully"));
    await session.step(117, "And the campaign \"BDD Submit\" should be saved with the name \"BDD Submit renamed\"", () => savedName(page, "BDD Submit", "BDD Submit renamed"));
    await session.step(118, "And \"Submit\" dialog should be visible", () => shouldBe(page, el("\"Submit\" dialog"), "visible"));
    await session.step(119, "And \"Name\" input in \"Submit\" dialog should have value \"BDD Submit renamed\"", () => shouldHaveValue(page, el("\"Name\" input in \"Submit\" dialog"), "BDD Submit renamed"));
    await session.step(120, "And \"Name\" input in \"Submit\" dialog should be read-only", () => shouldBe(page, el("\"Name\" input in \"Submit\" dialog"), "read-only"));
  });
  test("The pencil opens the file path, a folder that does not exist is refused, and the download icon saves the rows", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the Hit Design campaign \"BDD Submit\" with the stages \"Idea, Design\", the molecule field \"Reference compound\" set to \"c1ccc2[nH]ccc2c1\" and a table of:", () => hitDesignCampaignWithMolecule(page, "BDD Submit", "Idea, Design", "Reference compound", "c1ccc2[nH]ccc2c1", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["CCN","Idea","","45.08"]]);
    await session.step(25, "When user opens the campaign \"BDD Submit\" by its link", () => openByLink(page, "BDD Submit"));
    await session.step(26, "And user clicks on \"Submit...\" button in toolbar", () => clickOn(page, el("\"Submit...\" button in toolbar")));
    await session.step(27, "Then \"Submit\" dialog should be visible", () => shouldBe(page, el("\"Submit\" dialog"), "visible"));
    await session.step(123, "Given user watches downloads", () => watchDownloads(page));
    await session.step(124, "When user clicks on \"Edit file path\" icon in \"Submit\" dialog", () => clickOn(page, el("\"Edit file path\" icon in \"Submit\" dialog")));
    await session.step(125, "Then \"File Path\" input in \"Submit\" dialog should be editable", () => shouldBe(page, el("\"File Path\" input in \"Submit\" dialog"), "editable"));
    await session.step(126, "When user types \"System.AppData/HitTriage/BDD no such folder\" into \"File Path\" input in \"Submit\" dialog", () => typeInto(page, "System.AppData/HitTriage/BDD no such folder", el("\"File Path\" input in \"Submit\" dialog")));
    await session.step(127, "And user clicks on \"Save\" icon in \"Submit\" dialog", () => clickOn(page, el("\"Save\" icon in \"Submit\" dialog")));
    await session.step(128, "Then \"File Path\" input in \"Submit\" dialog should have the error \"Given folder does not exist\"", () => shouldHaveError(page, el("\"File Path\" input in \"Submit\" dialog"), "Given folder does not exist"));
    await session.step(129, "And an error balloon containing \"Given folder does not exist\" should have been shown", () => errorBalloonText(page, "Given folder does not exist"));
    await session.step(130, "When user clicks on \"Cancel\" icon in \"Submit\" dialog", () => clickOn(page, el("\"Cancel\" icon in \"Submit\" dialog")));
    await session.step(131, "Then \"File Path\" input in \"Submit\" dialog should have value \"System.AppData/HitTriage/Hit Design/campaigns\"", () => shouldHaveValue(page, el("\"File Path\" input in \"Submit\" dialog"), "System.AppData/HitTriage/Hit Design/campaigns"));
    await session.step(132, "And \"File Path\" input in \"Submit\" dialog should be read-only", () => shouldBe(page, el("\"File Path\" input in \"Submit\" dialog"), "read-only"));
    await session.step(133, "When user downloads a file through \"Download\" icon in \"Submit\" dialog", () => downloadThrough(page, el("\"Download\" icon in \"Submit\" dialog")));
    await session.step(134, "Then a file \"BDDBDDSUBMIT-1.csv\" should have been downloaded", () => fileDownloaded(page, "BDDBDDSUBMIT-1.csv"));
    await session.step(135, "And the downloaded file should contain \"CCN\"", () => downloadedContains(page, "CCN"));
  });
});
