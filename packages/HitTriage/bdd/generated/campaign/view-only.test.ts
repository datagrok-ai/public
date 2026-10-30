/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/campaign/view-only.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noLockHeld, openByLink, rememberUpdateLog, savedCell, savedVersion, updateLogAsRemembered, viewOnlyCampaign} from '../../bindings/campaign.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {doubleClickArea, errorBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A Hit Design campaign the user may only view", () => {
  const session = feature(test, "features/campaign/view-only.feature", import.meta.url);
  test("The campaign opens without the permissions icon, and an edit is refused and kept in the view", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the Hit Design campaign \"BDD View Only\", which only the sharing user may edit, with a table of:", () => viewOnlyCampaign(page, "BDD View Only", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"]]);
    await session.step(17, "When user opens the campaign \"BDD View Only\" by its link", () => openByLink(page, "BDD View Only"));
    await session.step(20, "Then \"Refresh\" icon in toolbar should be visible", () => shouldBe(page, el("\"Refresh\" icon in toolbar"), "visible"));
    await session.step(21, "And \"Submit...\" button in toolbar should be visible", () => shouldBe(page, el("\"Submit...\" button in toolbar"), "visible"));
    await session.step(22, "And \"Edit campaign permissions\" icon in toolbar should be absent", () => shouldBe(page, el("\"Edit campaign permissions\" icon in toolbar"), "absent"));
    await session.step(23, "When user double-clicks on the \"cell 1 of Stage\" area of grid", () => doubleClickArea(page, "cell 1 of Stage", el("grid")));
    await session.step(24, "And user selects \"Design\" in cell editor", () => selectIn(page, "Design", el("cell editor")));
    await session.step(25, "Then an error balloon containing \"You do not have permission to modify this campaign\" should have been shown", () => errorBalloonText(page, "You do not have permission to modify this campaign"));
    await session.step(26, "And the campaign \"BDD View Only\" should hold no lock on the server", () => noLockHeld(page, "BDD View Only"));
    await session.step(27, "And the campaign \"BDD View Only\" should be at version 1 on the server", () => savedVersion(page, "BDD View Only", 1));
    await session.step(28, "And the campaign \"BDD View Only\" should be saved with \"Idea\" in row 1 of the \"Stage\" column", () => savedCell(page, "BDD View Only", "Idea", 1, "Stage"));
    await session.step(29, "And the value of \"Stage\" column in row 1 should be \"Design\"", () => valueInRow(page, "Stage", 1, "Design"));
  });
  test("A refused edit logs no update, so no session is told of new changes", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the Hit Design campaign \"BDD View Only\", which only the sharing user may edit, with a table of:", () => viewOnlyCampaign(page, "BDD View Only", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"]]);
    await session.step(17, "When user opens the campaign \"BDD View Only\" by its link", () => openByLink(page, "BDD View Only"));
    await session.step(32, "When user remembers the update log of the campaign \"BDD View Only\"", () => rememberUpdateLog(page, "BDD View Only"));
    await session.step(33, "And user double-clicks on the \"cell 1 of Stage\" area of grid", () => doubleClickArea(page, "cell 1 of Stage", el("grid")));
    await session.step(34, "And user selects \"Design\" in cell editor", () => selectIn(page, "Design", el("cell editor")));
    await session.step(35, "Then an error balloon containing \"You do not have permission to modify this campaign\" should have been shown", () => errorBalloonText(page, "You do not have permission to modify this campaign"));
    await session.step(36, "And the campaign \"BDD View Only\" should hold no lock on the server", () => noLockHeld(page, "BDD View Only"));
    await session.step(37, "And the update log of the campaign \"BDD View Only\" should be as remembered", () => updateLogAsRemembered(page, "BDD View Only"));
  });
});
