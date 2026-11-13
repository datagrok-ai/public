/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/campaign/refresh-hint.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {hitDesignCampaign, openByLink, otherSessionSaves} from '../../bindings/campaign.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {elementHinted, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A Hit Design campaign view tells when another session saved the campaign", () => {
  const session = feature(test, "features/campaign/refresh-hint.feature", import.meta.url);
  test("Another session's save lights the Refresh icon, and Refresh brings its row in", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the Hit Design campaign \"BDD Refresh\" with the stages \"Idea, Design\" and a table of:", () => hitDesignCampaign(page, "BDD Refresh", "Idea, Design", [["Molecule","Stage","V-iD","MW"],["CCO","Design","","46.07"],["c1ccccc1","Design","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Design","","46.07"],["c1ccccc1","Design","","78.11"]]);
    await session.step(17, "When user opens the campaign \"BDD Refresh\" by its link", () => openByLink(page, "BDD Refresh"));
    await session.step(20, "When user hovers over \"Refresh\" icon in toolbar", () => hoverOver(page, el("\"Refresh\" icon in toolbar")));
    await session.step(21, "Then tooltip should have text \"Refresh Campaign\"", () => shouldHaveText(page, el("tooltip"), "Refresh Campaign"));
    await session.step(22, "When another session saves the campaign \"BDD Refresh\" with a new row of \"CCN\"", () => otherSessionSaves(page, "BDD Refresh", "CCN"));
    await session.step(23, "Then \"Refresh\" icon in toolbar should be hinted", () => elementHinted(page, el("\"Refresh\" icon in toolbar")));
    await session.step(24, "When user hovers over \"Refresh\" icon in toolbar", () => hoverOver(page, el("\"Refresh\" icon in toolbar")));
    await session.step(25, "Then tooltip should have text \"New changes available. Click to refresh campaign.\"", () => shouldHaveText(page, el("tooltip"), "New changes available. Click to refresh campaign."));
    await session.step(26, "And the table should have 2 rows", () => rowCount(page, 2));
    await session.step(27, "When user clicks on \"Refresh\" icon in toolbar", () => clickOn(page, el("\"Refresh\" icon in toolbar")));
    await session.step(28, "Then the table should have 3 rows", () => rowCount(page, 3));
    await session.step(29, "And the \"BDD Refresh\" view should be current", () => viewIsCurrent(page, "BDD Refresh"));
    await session.step(30, "And the value of \"Molecule\" column in row 3 should be \"CCN\"", () => valueInRow(page, "Molecule", 3, "CCN"));
    await session.step(31, "And the value of \"Stage\" column in row 3 should be \"Idea\"", () => valueInRow(page, "Stage", 3, "Idea"));
    await session.step(32, "When user hovers over \"Refresh\" icon in toolbar", () => hoverOver(page, el("\"Refresh\" icon in toolbar")));
    await session.step(33, "Then tooltip should have text \"Refresh Campaign\"", () => shouldHaveText(page, el("tooltip"), "Refresh Campaign"));
  });
});
