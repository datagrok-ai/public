/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/campaign/stage-choices.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {hitDesignCampaign, openByLink, savedCell, savedVersion} from '../../bindings/campaign.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey, selectIn, shouldBe, shouldHaveValue, shouldOffer} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {doubleClickArea, noBalloons, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Stage cell of a Hit Design campaign is a choice of the template's stages", () => {
  const session = feature(test, "features/campaign/stage-choices.feature", import.meta.url);
  test("The Stage cell offers the template's stages, and the stage picked is saved with no balloon", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the Hit Design campaign \"BDD Stage Choices\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Stage Choices", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Mkae","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Mkae","","78.11"]]);
    await session.step(17, "When user opens the campaign \"BDD Stage Choices\" by its link", () => openByLink(page, "BDD Stage Choices"));
    await session.step(18, "Then the table should have 2 rows", () => rowCount(page, 2));
    await session.step(21, "When user double-clicks on the \"cell 1 of Stage\" area of grid", () => doubleClickArea(page, "cell 1 of Stage", el("grid")));
    await session.step(22, "Then cell editor should offer \"Idea, Design, Make\"", () => shouldOffer(page, el("cell editor"), "Idea, Design, Make"));
    await session.step(23, "When user selects \"Make\" in cell editor", () => selectIn(page, "Make", el("cell editor")));
    await session.step(24, "Then the campaign \"BDD Stage Choices\" should be saved with \"Make\" in row 1 of the \"Stage\" column", () => savedCell(page, "BDD Stage Choices", "Make", 1, "Stage"));
    await session.step(25, "And the campaign \"BDD Stage Choices\" should be at version 2 on the server", () => savedVersion(page, "BDD Stage Choices", 2));
    await session.step(26, "And the value of \"Stage\" column in row 1 should be \"Make\"", () => valueInRow(page, "Stage", 1, "Make"));
    await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A stage the template does not have still shows, in the grid and in its dropdown, and Enter keeps it", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the Hit Design campaign \"BDD Stage Choices\" with the stages \"Idea, Design, Make\" and a table of:", () => hitDesignCampaign(page, "BDD Stage Choices", "Idea, Design, Make", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Mkae","","78.11"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"],["c1ccccc1","Mkae","","78.11"]]);
    await session.step(17, "When user opens the campaign \"BDD Stage Choices\" by its link", () => openByLink(page, "BDD Stage Choices"));
    await session.step(18, "Then the table should have 2 rows", () => rowCount(page, 2));
    await session.step(30, "Then the \"text of cell 2 of Stage\" reading of grid should be \"Mkae\"", () => readingReads(page, "text of cell 2 of Stage", el("grid"), "Mkae"));
    await session.step(31, "When user double-clicks on the \"cell 2 of Stage\" area of grid", () => doubleClickArea(page, "cell 2 of Stage", el("grid")));
    await session.step(32, "Then cell editor should offer \"Idea, Design, Make, Mkae\"", () => shouldOffer(page, el("cell editor"), "Idea, Design, Make, Mkae"));
    await session.step(33, "And cell editor should have the value \"Mkae\"", () => shouldHaveValue(page, el("cell editor"), "Mkae"));
    await session.step(34, "When user presses Enter", () => pressKey(page, "Enter"));
    await session.step(35, "Then cell editor should be absent", () => shouldBe(page, el("cell editor"), "absent"));
    await session.step(36, "And the campaign \"BDD Stage Choices\" should be at version 2 on the server", () => savedVersion(page, "BDD Stage Choices", 2));
    await session.step(37, "And the campaign \"BDD Stage Choices\" should be saved with \"Mkae\" in row 2 of the \"Stage\" column", () => savedCell(page, "BDD Stage Choices", "Mkae", 2, "Stage"));
    await session.step(38, "And the value of \"Stage\" column in row 2 should be \"Mkae\"", () => valueInRow(page, "Stage", 2, "Mkae"));
    await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
