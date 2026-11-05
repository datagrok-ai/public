/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/campaign/molecule-vid.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {enterMoleculeInSketcher, hitDesignCampaign, noLockHeld, openByLink, savedCellMatching, valueMatching} from '../../bindings/campaign.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, finishedUpdating, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentColumnIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShowsCurrentCell, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A molecule edited in a Hit Design campaign gets its V-iD", () => {
  const session = feature(test, "features/campaign/molecule-vid.feature", import.meta.url);
  test("A molecule typed into the sketcher gets a V-iD and the template's property, is saved, and its V-iD's pane lists the campaign", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(14, "And the Hit Design campaign \"BDD Vid\" with the stages \"Idea, Design\" and a table of:", () => hitDesignCampaign(page, "BDD Vid", "Idea, Design", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"]]);
    await session.step(17, "When user opens the campaign \"BDD Vid\" by its link", () => openByLink(page, "BDD Vid"));
    await session.step(20, "When user double-clicks on the \"cell 1 of Molecule\" area of grid", () => doubleClickArea(page, "cell 1 of Molecule", el("grid")));
    await session.step(21, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(22, "When user enters the molecule \"NCCCCCCCCCCCCCCCCCCCCCCCCO\" into sketcher dialog", () => enterMoleculeInSketcher(page, "NCCCCCCCCCCCCCCCCCCCCCCCCO"));
    await session.step(23, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(24, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(25, "And the campaign \"BDD Vid\" should be saved with a value matching \"^V[0-9]{6}$\" in row 1 of the \"V-iD\" column", () => savedCellMatching(page, "BDD Vid", "^V[0-9]{6}$", 1, "V-iD"));
    await session.step(26, "And the value of \"V-iD\" column in row 1 should match \"^V[0-9]{6}$\"", () => valueMatching(page, "V-iD", 1, "^V[0-9]{6}$"));
    await session.step(27, "And the value of \"MW\" column in row 1 should match \"^369[.]6\"", () => valueMatching(page, "MW", 1, "^369[.]6"));
    await session.step(28, "And the campaign \"BDD Vid\" should be saved with a value matching \"^369[.]6\" in row 1 of the \"MW\" column", () => savedCellMatching(page, "BDD Vid", "^369[.]6", 1, "MW"));
    await session.step(29, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(30, "And the campaign \"BDD Vid\" should hold no lock on the server", () => noLockHeld(page, "BDD Vid"));
    await session.step(31, "And grid should have finished updating", () => finishedUpdating(page, el("grid")));
    await session.step(32, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(33, "When user clicks on the \"cell 1 of V-iD\" area of grid", () => clickArea(page, "cell 1 of V-iD", el("grid")));
    await session.step(34, "Then the current column should be \"V-iD\"", () => currentColumnIs(page, "V-iD"));
    await session.step(35, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
    await session.step(36, "And \"Hit Design V-iD\" section should be visible", () => shouldBe(page, el("\"Hit Design V-iD\" section"), "visible"));
    await session.step(37, "When user expands \"Hit Design V-iD\" section", () => expand(page, el("\"Hit Design V-iD\" section")));
    await session.step(38, "Then \"Hit Design V-iD\" section should contain text \"Hit Design Campaigns\"", () => shouldContainText(page, el("\"Hit Design V-iD\" section"), "Hit Design Campaigns"));
    await session.step(39, "And \"Hit Design V-iD\" section should contain text \"BDDBDDVID-1\"", () => shouldContainText(page, el("\"Hit Design V-iD\" section"), "BDDBDDVID-1"));
    await session.step(40, "And \"Hit Design V-iD\" section should contain text \"{login}\"", () => shouldContainText(page, el("\"Hit Design V-iD\" section"), session.text("{login}")));
  });
});
