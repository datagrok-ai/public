/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/campaign/tracker-cards.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {cardLayoutOver, hitDesignCampaign, openByLink} from '../../bindings/campaign.js';
import {peptiHitCampaign} from '../../bindings/pepti-hit.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The progress tracker's cards are the campaign's own rows", () => {
  const session = feature(test, "features/campaign/tracker-cards.feature", import.meta.url);
  test("Next to a PeptiHit campaign, the cards of a kept card layout show this campaign's rows", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "Given the PeptiHit campaign \"BDD Cards Peptide\" with the stages \"Design\" and a table of:", () => peptiHitCampaign(page, "BDD Cards Peptide", "Design", [["Helm","Molecule","Stage","V-iD"],["PEPTIDE1{A.G}$$$$V2.0","CC(N)C(=O)NCC(=O)O","Design","V9900009"]]), [["Helm","Molecule","Stage","V-iD"],["PEPTIDE1{A.G}$$$$V2.0","CC(N)C(=O)NCC(=O)O","Design","V9900009"]]);
    await session.step(21, "And the Hit Design campaign \"BDD Cards\" with the stages \"Design, Make, Test\" and a table of:", () => hitDesignCampaign(page, "BDD Cards", "Design, Make, Test", [["Molecule","Stage","V-iD","MW"],["CCO","Design","V9900001","46.07"],["c1ccccc1","Make","V9900002","78.11"],["CCN","Design","V9900003","45.08"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Design","V9900001","46.07"],["c1ccccc1","Make","V9900002","78.11"],["CCN","Design","V9900003","45.08"]]);
    await session.step(26, "And the card layout of the campaign \"BDD Cards\" was made over a table named \"Table\"", () => cardLayoutOver(page, "BDD Cards", "Table"));
    await session.step(27, "When user opens the campaign \"BDD Cards Peptide\" by its link", () => openByLink(page, "BDD Cards Peptide"));
    await session.step(28, "And user opens the campaign \"BDD Cards\" by its link", () => openByLink(page, "BDD Cards"));
    await session.step(29, "And user clicks on \"Progress tracker\" button in toolbar", () => clickOn(page, el("\"Progress tracker\" button in toolbar")));
    await session.step(30, "Then \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "visible"));
    await session.step(31, "And the 'tiles in lane \"Design\"' reading of \"Tile Viewer\" viewer should be 2", () => readingIs(page, "tiles in lane \"Design\"", el("\"Tile Viewer\" viewer"), 2));
    await session.step(32, "And the \"V-iD of row 1\" reading of \"Tile Viewer\" viewer should be \"V9900001\"", () => readingReads(page, "V-iD of row 1", el("\"Tile Viewer\" viewer"), "V9900001"));
    await session.step(33, "And the \"V-iD of row 3\" reading of \"Tile Viewer\" viewer should be \"V9900003\"", () => readingReads(page, "V-iD of row 3", el("\"Tile Viewer\" viewer"), "V9900003"));
    await session.step(34, "And the \"lane of row 2\" reading of \"Tile Viewer\" viewer should be \"Make\"", () => readingReads(page, "lane of row 2", el("\"Tile Viewer\" viewer"), "Make"));
    await session.step(35, "And the \"V-iD of row 2\" reading of \"Tile Viewer\" viewer should be \"V9900002\"", () => readingReads(page, "V-iD of row 2", el("\"Tile Viewer\" viewer"), "V9900002"));
  });
  test("A row in a stage the template does not have is shown in a lane of that stage", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(38, "Given the Hit Design campaign \"BDD Odd Stage\" with the stages \"Design, Make, Test\" and a table of:", () => hitDesignCampaign(page, "BDD Odd Stage", "Design, Make, Test", [["Molecule","Stage","V-iD","MW"],["CCO","Design","","46.07"],["c1ccccc1","Mkae","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Design","","46.07"],["c1ccccc1","Mkae","","78.11"]]);
    await session.step(42, "When user opens the campaign \"BDD Odd Stage\" by its link", () => openByLink(page, "BDD Odd Stage"));
    await session.step(43, "And user clicks on \"Progress tracker\" button in toolbar", () => clickOn(page, el("\"Progress tracker\" button in toolbar")));
    await session.step(44, "Then \"Progress Tracker\" dialog should be visible", () => shouldBe(page, el("\"Progress Tracker\" dialog"), "visible"));
    await session.step(45, "And the \"lane names\" reading of \"Tile Viewer\" viewer should be \"Design, Make, Test, Mkae\"", () => readingReads(page, "lane names", el("\"Tile Viewer\" viewer"), "Design, Make, Test, Mkae"));
    await session.step(46, "And the \"lane of row 2\" reading of \"Tile Viewer\" viewer should be \"Mkae\"", () => readingReads(page, "lane of row 2", el("\"Tile Viewer\" viewer"), "Mkae"));
  });
});
