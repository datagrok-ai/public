/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/campaign/ribbon.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {hitDesignCampaign, openByLink, ribbonLastGroup} from '../../bindings/campaign.js';
import {notClipped, windowSize} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, toolboxPaneShown, urlShouldContain} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {menuLists, noBalloons, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Hit Design campaign view's ribbon and crumbs", () => {
  const session = feature(test, "features/campaign/ribbon.feature", import.meta.url);
  test("The ribbon ends with the package's actions, and the crumbs name the app and the campaign", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the Hit Design campaign \"BDD Ribbon\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Ribbon", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]);
    await session.step(21, "When user opens the campaign \"BDD Ribbon\" by its link", () => openByLink(page, "BDD Ribbon"));
    await session.step(24, "Then the last group of the ribbon should list \"Refresh, Add new row, Merge another table into the campaign, Calculate additional properties, Edit campaign permissions, Apply template layout (Loaded from mounted file storage), Progress tracker, Submit...\"", () => ribbonLastGroup(page, "Refresh, Add new row, Merge another table into the campaign, Calculate additional properties, Edit campaign permissions, Apply template layout (Loaded from mounted file storage), Progress tracker, Submit..."));
    await session.step(25, "And \"Save\" button in toolbar should be hidden", () => shouldBe(page, el("\"Save\" button in toolbar"), "hidden"));
    await session.step(26, "And \"Hit Design\" breadcrumb should be visible", () => shouldBe(page, el("\"Hit Design\" breadcrumb"), "visible"));
    await session.step(27, "And \"BDD Ribbon\" breadcrumb should be visible", () => shouldBe(page, el("\"BDD Ribbon\" breadcrumb"), "visible"));
    await session.step(28, "And the page address should contain \"campaignId=BDDBDDRIBBON-1\"", () => urlShouldContain(page, "campaignId=BDDBDDRIBBON-1"));
    await session.step(29, "When user hovers over \"Refresh\" icon in toolbar", () => hoverOver(page, el("\"Refresh\" icon in toolbar")));
    await session.step(30, "Then tooltip should have text \"Refresh Campaign\"", () => shouldHaveText(page, el("tooltip"), "Refresh Campaign"));
    await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Apply template layout loads the template's layout file into the view", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the Hit Design campaign \"BDD Ribbon\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Ribbon", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]);
    await session.step(21, "When user opens the campaign \"BDD Ribbon\" by its link", () => openByLink(page, "BDD Ribbon"));
    await session.step(34, "Then the open tableview should have 0 Bar chart viewers", () => viewerCount(page, 0, "Bar chart"));
    await session.step(35, "When user clicks on \"Apply template layout (Loaded from mounted file storage)\" icon in toolbar", () => clickOn(page, el("\"Apply template layout (Loaded from mounted file storage)\" icon in toolbar")));
    await session.step(36, "Then the open tableview should have 4 Bar chart viewers", () => viewerCount(page, 4, "Bar chart"));
    await session.step(37, "And the table should have 2 rows", () => rowCount(page, 2));
  });
  test("At 1366 px the primary buttons stay in view and the package's icons are one menu away", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the Hit Design campaign \"BDD Ribbon\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Ribbon", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Design","","78.11"]]);
    await session.step(21, "When user opens the campaign \"BDD Ribbon\" by its link", () => openByLink(page, "BDD Ribbon"));
    await session.step(40, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(41, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(42, "When the window is 1366 by 768 px", () => windowSize(page, 1366, 768));
    await session.step(43, "Then \"Submit...\" button in toolbar should not be clipped", () => notClipped(page, el("\"Submit...\" button in toolbar")));
    await session.step(44, "And \"Progress tracker\" button in toolbar should not be clipped", () => notClipped(page, el("\"Progress tracker\" button in toolbar")));
    await session.step(45, "And \"Merge another table into the campaign\" icon in toolbar should be hidden", () => shouldBe(page, el("\"Merge another table into the campaign\" icon in toolbar"), "hidden"));
    await session.step(46, "When user clicks on \"More actions\" icon in toolbar", () => clickOn(page, el("\"More actions\" icon in toolbar")));
    await session.step(47, "Then the open menu should list \"Refresh Campaign\"", () => menuLists(page, "Refresh Campaign"));
    await session.step(48, "And the open menu should list \"Merge another table into the campaign...\"", () => menuLists(page, "Merge another table into the campaign..."));
    await session.step(49, "And the open menu should list \"Edit campaign permissions...\"", () => menuLists(page, "Edit campaign permissions..."));
  });
});
