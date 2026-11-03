/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/campaign/grid.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {hitDesignCampaign, openByLink, savedCellMatching, valueMatching} from '../../bindings/campaign.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeContextMenu, menuDoesNotList, menuLists, noBalloons, pickFromAreaContextMenu, rightClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Rows and cells of a Hit Design campaign", () => {
  const session = feature(test, "features/campaign/grid.feature", import.meta.url);
  test("Add new row gives the new row the first stage and opens the sketcher on its empty molecule", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Hit Design campaign \"BDD Grid\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Grid", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Design","","46.07"],["c1ccccc1","Design","",""]]), [["Molecule","Stage","V-iD","MW"],["CCO","Design","","46.07"],["c1ccccc1","Design","",""]]);
    await session.step(20, "When user opens the campaign \"BDD Grid\" by its link", () => openByLink(page, "BDD Grid"));
    await session.step(21, "Then the table should have 2 rows", () => rowCount(page, 2));
    await session.step(24, "When user clicks on \"Add new row\" icon in toolbar", () => clickOn(page, el("\"Add new row\" icon in toolbar")));
    await session.step(25, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(26, "And the table should have 3 rows", () => rowCount(page, 3));
    await session.step(27, "And the value of \"Stage\" column in row 3 should be \"Idea\"", () => valueInRow(page, "Stage", 3, "Idea"));
    await session.step(28, "And the value of \"Stage\" column in row 1 should be \"Design\"", () => valueInRow(page, "Stage", 1, "Design"));
    await session.step(29, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
    await session.step(30, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(31, "And the table should have 3 rows", () => rowCount(page, 3));
  });
  test("A Molecule cell's menu adds the package's row actions, and another column's menu does not", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Hit Design campaign \"BDD Grid\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Grid", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Design","","46.07"],["c1ccccc1","Design","",""]]), [["Molecule","Stage","V-iD","MW"],["CCO","Design","","46.07"],["c1ccccc1","Design","",""]]);
    await session.step(20, "When user opens the campaign \"BDD Grid\" by its link", () => openByLink(page, "BDD Grid"));
    await session.step(21, "Then the table should have 2 rows", () => rowCount(page, 2));
    await session.step(34, "When user right-clicks on the \"cell 1 of Molecule\" area of grid", () => rightClickArea(page, "cell 1 of Molecule", el("grid")));
    await session.step(35, "Then the open menu should list \"Add new row\"", () => menuLists(page, "Add new row"));
    await session.step(36, "And the open menu should list \"Duplicate molecule\"", () => menuLists(page, "Duplicate molecule"));
    await session.step(37, "And the open menu should list \"Re-Run Calculations\"", () => menuLists(page, "Re-Run Calculations"));
    await session.step(38, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(39, "And user right-clicks on the \"cell 1 of Stage\" area of grid", () => rightClickArea(page, "cell 1 of Stage", el("grid")));
    await session.step(40, "Then the open menu should list \"Sort...\"", () => menuLists(page, "Sort..."));
    await session.step(41, "And the open menu should not list \"Duplicate molecule\"", () => menuDoesNotList(page, "Duplicate molecule"));
    await session.step(42, "And the open menu should not list \"Re-Run Calculations\"", () => menuDoesNotList(page, "Re-Run Calculations"));
  });
  test("Duplicate molecule appends the cell's molecule as a new row in the first stage", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Hit Design campaign \"BDD Grid\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Grid", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Design","","46.07"],["c1ccccc1","Design","",""]]), [["Molecule","Stage","V-iD","MW"],["CCO","Design","","46.07"],["c1ccccc1","Design","",""]]);
    await session.step(20, "When user opens the campaign \"BDD Grid\" by its link", () => openByLink(page, "BDD Grid"));
    await session.step(21, "Then the table should have 2 rows", () => rowCount(page, 2));
    await session.step(45, "When user picks \"Duplicate molecule\" from the context menu of the \"cell 1 of Molecule\" area of grid", () => pickFromAreaContextMenu(page, "Duplicate molecule", "cell 1 of Molecule", el("grid")));
    await session.step(46, "Then the table should have 3 rows", () => rowCount(page, 3));
    await session.step(47, "And the value of \"Molecule\" column in row 3 should be \"CCO\"", () => valueInRow(page, "Molecule", 3, "CCO"));
    await session.step(48, "And the value of \"Stage\" column in row 3 should be \"Idea\"", () => valueInRow(page, "Stage", 3, "Idea"));
  });
  test("Re-Run Calculations computes the template's property for that row alone and saves it", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Hit Design campaign \"BDD Grid\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Grid", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Design","","46.07"],["c1ccccc1","Design","",""]]), [["Molecule","Stage","V-iD","MW"],["CCO","Design","","46.07"],["c1ccccc1","Design","",""]]);
    await session.step(20, "When user opens the campaign \"BDD Grid\" by its link", () => openByLink(page, "BDD Grid"));
    await session.step(21, "Then the table should have 2 rows", () => rowCount(page, 2));
    await session.step(51, "Then the value of \"MW\" column in row 2 should be \"\"", () => valueInRow(page, "MW", 2, ""));
    await session.step(52, "When user picks \"Re-Run Calculations\" from the context menu of the \"cell 2 of Molecule\" area of grid", () => pickFromAreaContextMenu(page, "Re-Run Calculations", "cell 2 of Molecule", el("grid")));
    await session.step(53, "Then the campaign \"BDD Grid\" should be saved with a value matching \"^78[.]1\" in row 2 of the \"MW\" column", () => savedCellMatching(page, "BDD Grid", "^78[.]1", 2, "MW"));
    await session.step(54, "And the value of \"MW\" column in row 2 should match \"^78[.]1\"", () => valueMatching(page, "MW", 2, "^78[.]1"));
    await session.step(55, "And the value of \"MW\" column in row 1 should match \"^46[.]0[67]\"", () => valueMatching(page, "MW", 1, "^46[.]0[67]"));
    await session.step(56, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
