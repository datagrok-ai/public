/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/campaign/merge.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {chooseProjectFile, enterUserFilePath, hitDesignCampaign, holdsUserFilePath, openByLink, savedCell, savedRows, savedToMergeNoFile, savedToMergeUserFile, userFileHolds, valueMatching} from '../../bindings/campaign.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, doubleClickOn, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {infoBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Merging a table from the user's files into a Hit Design campaign", () => {
  const session = feature(test, "features/campaign/merge.feature", import.meta.url);
  test("OK stays disabled, saying why, until the dialog has what the merge needs", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the Hit Design campaign \"BDD Merge\" with the stages \"Idea, Design\" and a table of:", () => hitDesignCampaign(page, "BDD Merge", "Idea, Design", [["Molecule","Stage","V-iD","MW"],["CCO","Design","V900001","46.07"],["c1ccccc1","Design","V900002","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Design","V900001","46.07"],["c1ccccc1","Design","V900002","78.11"]]);
    await session.step(23, "And the user's file \"bdd-merge-by-smiles.csv\" holds the table:", () => userFileHolds(page, "bdd-merge-by-smiles.csv", [["Molecule","Score"],["OCC","5"],["NCCCCCCCCCCCCCCCCCCCCCCCCO","7"]]), [["Molecule","Score"],["OCC","5"],["NCCCCCCCCCCCCCCCCCCCCCCCCO","7"]]);
    await session.step(27, "And the user's file \"bdd-merge-by-vid.csv\" holds the table:", () => userFileHolds(page, "bdd-merge-by-vid.csv", [["V-iD","Priority"],["V900002","High"]]), [["V-iD","Priority"],["V900002","High"]]);
    await session.step(30, "When user opens the campaign \"BDD Merge\" by its link", () => openByLink(page, "BDD Merge"));
    await session.step(31, "And user clicks on \"Merge another table into the campaign\" icon in toolbar", () => clickOn(page, el("\"Merge another table into the campaign\" icon in toolbar")));
    await session.step(32, "Then \"Merge with table\" dialog should be visible", () => shouldBe(page, el("\"Merge with table\" dialog"), "visible"));
    await session.step(35, "Then OK button in \"Merge with table\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Merge with table\" dialog"), "disabled"));
    await session.step(36, "And \"Merge with table\" dialog should contain text \"Pick a source file to merge from.\"", () => shouldContainText(page, el("\"Merge with table\" dialog"), "Pick a source file to merge from."));
    await session.step(37, "And CANCEL button in \"Merge with table\" dialog should be enabled", () => shouldBe(page, el("CANCEL button in \"Merge with table\" dialog"), "enabled"));
    await session.step(38, "When user enters the path of the user's file \"bdd-merge-by-smiles.csv\" into \"Source file\" input in \"Merge with table\" dialog", () => enterUserFilePath(page, "bdd-merge-by-smiles.csv", el("\"Source file\" input in \"Merge with table\" dialog")));
    await session.step(39, "Then \"Merge with table\" dialog should contain text \"Pick the V-iD column from the incoming table.\"", () => shouldContainText(page, el("\"Merge with table\" dialog"), "Pick the V-iD column from the incoming table."));
    await session.step(40, "And OK button in \"Merge with table\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Merge with table\" dialog"), "disabled"));
    await session.step(41, "When user selects \"By canonical SMILES\" in \"Merge by\" input in \"Merge with table\" dialog", () => selectIn(page, "By canonical SMILES", el("\"Merge by\" input in \"Merge with table\" dialog")));
    await session.step(42, "Then OK button in \"Merge with table\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Merge with table\" dialog"), "enabled"));
    await session.step(43, "When user selects \"(none)\" in \"Molecule column (incoming)\" input in \"Merge with table\" dialog", () => selectIn(page, "(none)", el("\"Molecule column (incoming)\" input in \"Merge with table\" dialog")));
    await session.step(44, "Then \"Merge with table\" dialog should contain text \"Pick the molecule column from the incoming table.\"", () => shouldContainText(page, el("\"Merge with table\" dialog"), "Pick the molecule column from the incoming table."));
    await session.step(45, "And OK button in \"Merge with table\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Merge with table\" dialog"), "disabled"));
    await session.step(46, "When user clicks on CANCEL button in \"Merge with table\" dialog", () => clickOn(page, el("CANCEL button in \"Merge with table\" dialog")));
    await session.step(47, "Then \"Merge with table\" dialog should be absent", () => shouldBe(page, el("\"Merge with table\" dialog"), "absent"));
    await session.step(48, "And the table should have 2 rows", () => rowCount(page, 2));
  });
  test("By canonical SMILES, the matching row takes the file's values and the new molecule comes in with a stage, a V-iD and MW", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the Hit Design campaign \"BDD Merge\" with the stages \"Idea, Design\" and a table of:", () => hitDesignCampaign(page, "BDD Merge", "Idea, Design", [["Molecule","Stage","V-iD","MW"],["CCO","Design","V900001","46.07"],["c1ccccc1","Design","V900002","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Design","V900001","46.07"],["c1ccccc1","Design","V900002","78.11"]]);
    await session.step(23, "And the user's file \"bdd-merge-by-smiles.csv\" holds the table:", () => userFileHolds(page, "bdd-merge-by-smiles.csv", [["Molecule","Score"],["OCC","5"],["NCCCCCCCCCCCCCCCCCCCCCCCCO","7"]]), [["Molecule","Score"],["OCC","5"],["NCCCCCCCCCCCCCCCCCCCCCCCCO","7"]]);
    await session.step(27, "And the user's file \"bdd-merge-by-vid.csv\" holds the table:", () => userFileHolds(page, "bdd-merge-by-vid.csv", [["V-iD","Priority"],["V900002","High"]]), [["V-iD","Priority"],["V900002","High"]]);
    await session.step(30, "When user opens the campaign \"BDD Merge\" by its link", () => openByLink(page, "BDD Merge"));
    await session.step(31, "And user clicks on \"Merge another table into the campaign\" icon in toolbar", () => clickOn(page, el("\"Merge another table into the campaign\" icon in toolbar")));
    await session.step(32, "Then \"Merge with table\" dialog should be visible", () => shouldBe(page, el("\"Merge with table\" dialog"), "visible"));
    await session.step(51, "When user enters the path of the user's file \"bdd-merge-by-smiles.csv\" into \"Source file\" input in \"Merge with table\" dialog", () => enterUserFilePath(page, "bdd-merge-by-smiles.csv", el("\"Source file\" input in \"Merge with table\" dialog")));
    await session.step(52, "And user selects \"By canonical SMILES\" in \"Merge by\" input in \"Merge with table\" dialog", () => selectIn(page, "By canonical SMILES", el("\"Merge by\" input in \"Merge with table\" dialog")));
    await session.step(53, "And user clicks on OK button in \"Merge with table\" dialog", () => clickOn(page, el("OK button in \"Merge with table\" dialog")));
    await session.step(54, "Then an info balloon containing \"Table merged into campaign\" should have been shown", () => infoBalloonText(page, "Table merged into campaign"));
    await session.step(55, "And the campaign \"BDD Merge\" should be saved with 3 rows", () => savedRows(page, "BDD Merge", 3));
    await session.step(56, "And the table should have 3 rows", () => rowCount(page, 3));
    await session.step(57, "And the value of \"Score\" column in row 1 should be \"5\"", () => valueInRow(page, "Score", 1, "5"));
    await session.step(58, "And the value of \"Score\" column in row 2 should be \"\"", () => valueInRow(page, "Score", 2, ""));
    await session.step(59, "And the value of \"Molecule\" column in row 3 should be \"NCCCCCCCCCCCCCCCCCCCCCCCCO\"", () => valueInRow(page, "Molecule", 3, "NCCCCCCCCCCCCCCCCCCCCCCCCO"));
    await session.step(60, "And the value of \"Stage\" column in row 3 should be \"Idea\"", () => valueInRow(page, "Stage", 3, "Idea"));
    await session.step(61, "And the value of \"V-iD\" column in row 3 should match \"^V[0-9]{6}$\"", () => valueMatching(page, "V-iD", 3, "^V[0-9]{6}$"));
    await session.step(62, "And the value of \"MW\" column in row 3 should match \"^369[.]6\"", () => valueMatching(page, "MW", 3, "^369[.]6"));
    await session.step(63, "And the campaign \"BDD Merge\" should be saved with \"5\" in row 1 of the \"Score\" column", () => savedCell(page, "BDD Merge", "5", 1, "Score"));
  });
  test("By V-iD, the row of that V-iD takes the file's values and no row is added", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the Hit Design campaign \"BDD Merge\" with the stages \"Idea, Design\" and a table of:", () => hitDesignCampaign(page, "BDD Merge", "Idea, Design", [["Molecule","Stage","V-iD","MW"],["CCO","Design","V900001","46.07"],["c1ccccc1","Design","V900002","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Design","V900001","46.07"],["c1ccccc1","Design","V900002","78.11"]]);
    await session.step(23, "And the user's file \"bdd-merge-by-smiles.csv\" holds the table:", () => userFileHolds(page, "bdd-merge-by-smiles.csv", [["Molecule","Score"],["OCC","5"],["NCCCCCCCCCCCCCCCCCCCCCCCCO","7"]]), [["Molecule","Score"],["OCC","5"],["NCCCCCCCCCCCCCCCCCCCCCCCCO","7"]]);
    await session.step(27, "And the user's file \"bdd-merge-by-vid.csv\" holds the table:", () => userFileHolds(page, "bdd-merge-by-vid.csv", [["V-iD","Priority"],["V900002","High"]]), [["V-iD","Priority"],["V900002","High"]]);
    await session.step(30, "When user opens the campaign \"BDD Merge\" by its link", () => openByLink(page, "BDD Merge"));
    await session.step(31, "And user clicks on \"Merge another table into the campaign\" icon in toolbar", () => clickOn(page, el("\"Merge another table into the campaign\" icon in toolbar")));
    await session.step(32, "Then \"Merge with table\" dialog should be visible", () => shouldBe(page, el("\"Merge with table\" dialog"), "visible"));
    await session.step(66, "When user enters the path of the user's file \"bdd-merge-by-vid.csv\" into \"Source file\" input in \"Merge with table\" dialog", () => enterUserFilePath(page, "bdd-merge-by-vid.csv", el("\"Source file\" input in \"Merge with table\" dialog")));
    await session.step(67, "Then \"V-iD column (incoming)\" input in \"Merge with table\" dialog should have value \"V-iD\"", () => shouldHaveValue(page, el("\"V-iD column (incoming)\" input in \"Merge with table\" dialog"), "V-iD"));
    await session.step(68, "And OK button in \"Merge with table\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Merge with table\" dialog"), "enabled"));
    await session.step(69, "When user clicks on OK button in \"Merge with table\" dialog", () => clickOn(page, el("OK button in \"Merge with table\" dialog")));
    await session.step(70, "Then an info balloon containing \"Table merged into campaign\" should have been shown", () => infoBalloonText(page, "Table merged into campaign"));
    await session.step(71, "And the campaign \"BDD Merge\" should be saved with \"High\" in row 2 of the \"Priority\" column", () => savedCell(page, "BDD Merge", "High", 2, "Priority"));
    await session.step(72, "And the table should have 2 rows", () => rowCount(page, 2));
    await session.step(73, "And the value of \"Priority\" column in row 2 should be \"High\"", () => valueInRow(page, "Priority", 2, "High"));
    await session.step(74, "And the value of \"Priority\" column in row 1 should be \"\"", () => valueInRow(page, "Priority", 1, ""));
  });
  test("A file picked from Files is merged and remembered with the campaign, and the dialog opens on it again", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the Hit Design campaign \"BDD Merge\" with the stages \"Idea, Design\" and a table of:", () => hitDesignCampaign(page, "BDD Merge", "Idea, Design", [["Molecule","Stage","V-iD","MW"],["CCO","Design","V900001","46.07"],["c1ccccc1","Design","V900002","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Design","V900001","46.07"],["c1ccccc1","Design","V900002","78.11"]]);
    await session.step(23, "And the user's file \"bdd-merge-by-smiles.csv\" holds the table:", () => userFileHolds(page, "bdd-merge-by-smiles.csv", [["Molecule","Score"],["OCC","5"],["NCCCCCCCCCCCCCCCCCCCCCCCCO","7"]]), [["Molecule","Score"],["OCC","5"],["NCCCCCCCCCCCCCCCCCCCCCCCCO","7"]]);
    await session.step(27, "And the user's file \"bdd-merge-by-vid.csv\" holds the table:", () => userFileHolds(page, "bdd-merge-by-vid.csv", [["V-iD","Priority"],["V900002","High"]]), [["V-iD","Priority"],["V900002","High"]]);
    await session.step(30, "When user opens the campaign \"BDD Merge\" by its link", () => openByLink(page, "BDD Merge"));
    await session.step(31, "And user clicks on \"Merge another table into the campaign\" icon in toolbar", () => clickOn(page, el("\"Merge another table into the campaign\" icon in toolbar")));
    await session.step(32, "Then \"Merge with table\" dialog should be visible", () => shouldBe(page, el("\"Merge with table\" dialog"), "visible"));
    await session.step(79, "When user clicks on \"Add file from Files\" icon in \"Source file\" input in \"Merge with table\" dialog", () => clickOn(page, el("\"Add file from Files\" icon in \"Source file\" input in \"Merge with table\" dialog")));
    await session.step(80, "Then \"Select a file\" dialog should be visible", () => shouldBe(page, el("\"Select a file\" dialog"), "visible"));
    await session.step(81, "Given \"My files\" tree node in \"Select a file\" dialog is expanded", () => isExpanded(page, el("\"My files\" tree node in \"Select a file\" dialog")));
    await session.step(82, "When user double-clicks on \"bdd-merge-by-vid.csv\" tree node in \"Select a file\" dialog", () => doubleClickOn(page, el("\"bdd-merge-by-vid.csv\" tree node in \"Select a file\" dialog")));
    await session.step(83, "Then \"Select a file\" dialog should be absent", () => shouldBe(page, el("\"Select a file\" dialog"), "absent"));
    await session.step(84, "And \"Source file\" input in \"Merge with table\" dialog should hold the path of the user's file \"bdd-merge-by-vid.csv\"", () => holdsUserFilePath(page, el("\"Source file\" input in \"Merge with table\" dialog"), "bdd-merge-by-vid.csv"));
    await session.step(85, "And \"V-iD column (incoming)\" input in \"Merge with table\" dialog should have value \"V-iD\"", () => shouldHaveValue(page, el("\"V-iD column (incoming)\" input in \"Merge with table\" dialog"), "V-iD"));
    await session.step(86, "When user checks \"Auto-merge on campaign open\" input in \"Merge with table\" dialog", () => check(page, el("\"Auto-merge on campaign open\" input in \"Merge with table\" dialog")));
    await session.step(87, "And user clicks on OK button in \"Merge with table\" dialog", () => clickOn(page, el("OK button in \"Merge with table\" dialog")));
    await session.step(88, "Then an info balloon containing \"Table merged into campaign\" should have been shown", () => infoBalloonText(page, "Table merged into campaign"));
    await session.step(89, "And the campaign \"BDD Merge\" should be saved with \"High\" in row 2 of the \"Priority\" column", () => savedCell(page, "BDD Merge", "High", 2, "Priority"));
    await session.step(90, "And the campaign \"BDD Merge\" should be saved to merge the user's file \"bdd-merge-by-vid.csv\" on every open", () => savedToMergeUserFile(page, "BDD Merge", "bdd-merge-by-vid.csv"));
    await session.step(91, "When user clicks on \"Merge another table into the campaign\" icon in toolbar", () => clickOn(page, el("\"Merge another table into the campaign\" icon in toolbar")));
    await session.step(92, "Then \"Source file\" input in \"Merge with table\" dialog should hold the path of the user's file \"bdd-merge-by-vid.csv\"", () => holdsUserFilePath(page, el("\"Source file\" input in \"Merge with table\" dialog"), "bdd-merge-by-vid.csv"));
    await session.step(93, "And \"V-iD column (incoming)\" input in \"Merge with table\" dialog should have value \"V-iD\"", () => shouldHaveValue(page, el("\"V-iD column (incoming)\" input in \"Merge with table\" dialog"), "V-iD"));
    await session.step(94, "And \"Auto-merge on campaign open\" input in \"Merge with table\" dialog should be checked", () => shouldBe(page, el("\"Auto-merge on campaign open\" input in \"Merge with table\" dialog"), "checked"));
  });
  test("A file uploaded from the computer is merged once, with Auto-merge on campaign open off", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the Hit Design campaign \"BDD Merge\" with the stages \"Idea, Design\" and a table of:", () => hitDesignCampaign(page, "BDD Merge", "Idea, Design", [["Molecule","Stage","V-iD","MW"],["CCO","Design","V900001","46.07"],["c1ccccc1","Design","V900002","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Design","V900001","46.07"],["c1ccccc1","Design","V900002","78.11"]]);
    await session.step(23, "And the user's file \"bdd-merge-by-smiles.csv\" holds the table:", () => userFileHolds(page, "bdd-merge-by-smiles.csv", [["Molecule","Score"],["OCC","5"],["NCCCCCCCCCCCCCCCCCCCCCCCCO","7"]]), [["Molecule","Score"],["OCC","5"],["NCCCCCCCCCCCCCCCCCCCCCCCCO","7"]]);
    await session.step(27, "And the user's file \"bdd-merge-by-vid.csv\" holds the table:", () => userFileHolds(page, "bdd-merge-by-vid.csv", [["V-iD","Priority"],["V900002","High"]]), [["V-iD","Priority"],["V900002","High"]]);
    await session.step(30, "When user opens the campaign \"BDD Merge\" by its link", () => openByLink(page, "BDD Merge"));
    await session.step(31, "And user clicks on \"Merge another table into the campaign\" icon in toolbar", () => clickOn(page, el("\"Merge another table into the campaign\" icon in toolbar")));
    await session.step(32, "Then \"Merge with table\" dialog should be visible", () => shouldBe(page, el("\"Merge with table\" dialog"), "visible"));
    await session.step(98, "When user clicks on \"Open file\" icon in \"Source file\" input in \"Merge with table\" dialog and chooses the \"fixtures/merge-local.csv\" file of the project", () => chooseProjectFile(page, el("\"Open file\" icon in \"Source file\" input in \"Merge with table\" dialog"), "fixtures/merge-local.csv"));
    await session.step(99, "Then \"V-iD column (incoming)\" input in \"Merge with table\" dialog should have value \"V-iD\"", () => shouldHaveValue(page, el("\"V-iD column (incoming)\" input in \"Merge with table\" dialog"), "V-iD"));
    await session.step(100, "And \"Auto-merge on campaign open\" input in \"Merge with table\" dialog should be disabled", () => shouldBe(page, el("\"Auto-merge on campaign open\" input in \"Merge with table\" dialog"), "disabled"));
    await session.step(101, "When user clicks on OK button in \"Merge with table\" dialog", () => clickOn(page, el("OK button in \"Merge with table\" dialog")));
    await session.step(102, "Then an info balloon containing \"Table merged into campaign\" should have been shown", () => infoBalloonText(page, "Table merged into campaign"));
    await session.step(103, "And the campaign \"BDD Merge\" should be saved with \"Urgent\" in row 1 of the \"Priority\" column", () => savedCell(page, "BDD Merge", "Urgent", 1, "Priority"));
    await session.step(104, "And the campaign \"BDD Merge\" should be saved to merge no file on open", () => savedToMergeNoFile(page, "BDD Merge"));
    await session.step(105, "And the value of \"Priority\" column in row 1 should be \"Urgent\"", () => valueInRow(page, "Priority", 1, "Urgent"));
  });
});
