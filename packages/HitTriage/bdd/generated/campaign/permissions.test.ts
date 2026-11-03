/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/campaign/permissions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {grantSharingUserEdit, groupPicture, hitDesignCampaign, namesSharingUser, openByLink, sharingUserMayEdit, sharingUserPicture} from '../../bindings/campaign.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {infoBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Editing a Hit Design campaign's permissions from its ribbon", () => {
  const session = feature(test, "features/campaign/permissions.feature", import.meta.url);
  test("A user given edit access is saved with the campaign and shown when the dialog opens again", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the Hit Design campaign \"BDD Permissions\" with the stages \"Idea\" and a table of:", () => hitDesignCampaign(page, "BDD Permissions", "Idea", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"]]);
    await session.step(18, "When user opens the campaign \"BDD Permissions\" by its link", () => openByLink(page, "BDD Permissions"));
    await session.step(21, "When user clicks on \"Edit campaign permissions\" icon in toolbar", () => clickOn(page, el("\"Edit campaign permissions\" icon in toolbar")));
    await session.step(22, "Then \"Edit Permissions\" dialog should be visible", () => shouldBe(page, el("\"Edit Permissions\" dialog"), "visible"));
    await session.step(23, "And \"Edit Permissions\" dialog should contain text \"Can View\"", () => shouldContainText(page, el("\"Edit Permissions\" dialog"), "Can View"));
    await session.step(24, "And \"Edit Permissions\" dialog should contain text \"Can Edit\"", () => shouldContainText(page, el("\"Edit Permissions\" dialog"), "Can Edit"));
    await session.step(25, "And \"Edit Permissions\" dialog should contain text \"All users\"", () => shouldContainText(page, el("\"Edit Permissions\" dialog"), "All users"));
    await session.step(26, "And CANCEL button in \"Edit Permissions\" dialog should be enabled", () => shouldBe(page, el("CANCEL button in \"Edit Permissions\" dialog"), "enabled"));
    await session.step(27, "When user lets the sharing user edit in \"Edit Permissions\" dialog", () => grantSharingUserEdit(page, el("\"Edit Permissions\" dialog")));
    await session.step(28, "And user clicks on OK button in \"Edit Permissions\" dialog", () => clickOn(page, el("OK button in \"Edit Permissions\" dialog")));
    await session.step(29, "Then \"Edit Permissions\" dialog should be absent", () => shouldBe(page, el("\"Edit Permissions\" dialog"), "absent"));
    await session.step(30, "And an info balloon containing \"Campaign saved successfully\" should have been shown", () => infoBalloonText(page, "Campaign saved successfully"));
    await session.step(31, "And the campaign \"BDD Permissions\" should let the sharing user edit it on the server", () => sharingUserMayEdit(page, "BDD Permissions"));
    await session.step(32, "When user clicks on \"Edit campaign permissions\" icon in toolbar", () => clickOn(page, el("\"Edit campaign permissions\" icon in toolbar")));
    await session.step(33, "Then \"Edit Permissions\" dialog should be visible", () => shouldBe(page, el("\"Edit Permissions\" dialog"), "visible"));
    await session.step(34, "And \"Edit Permissions\" dialog should name the sharing user", () => namesSharingUser(page, el("\"Edit Permissions\" dialog")));
    await session.step(35, "And the sharing user's row of \"Edit Permissions\" dialog should show the user's picture", () => sharingUserPicture(page, el("\"Edit Permissions\" dialog")));
    await session.step(36, "When user clicks on CANCEL button in \"Edit Permissions\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Permissions\" dialog")));
    await session.step(37, "Then \"Edit Permissions\" dialog should be absent", () => shouldBe(page, el("\"Edit Permissions\" dialog"), "absent"));
  });
  test("A group row shows the group icon, a user added from the picker the user's picture", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the Hit Design campaign \"BDD Permissions\" with the stages \"Idea\" and a table of:", () => hitDesignCampaign(page, "BDD Permissions", "Idea", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"]]);
    await session.step(18, "When user opens the campaign \"BDD Permissions\" by its link", () => openByLink(page, "BDD Permissions"));
    await session.step(40, "When user clicks on \"Edit campaign permissions\" icon in toolbar", () => clickOn(page, el("\"Edit campaign permissions\" icon in toolbar")));
    await session.step(41, "Then \"Edit Permissions\" dialog should be visible", () => shouldBe(page, el("\"Edit Permissions\" dialog"), "visible"));
    await session.step(42, "And the \"All users\" row of \"Edit Permissions\" dialog should show a group icon", () => groupPicture(page, "All users", el("\"Edit Permissions\" dialog")));
    await session.step(43, "When user lets the sharing user edit in \"Edit Permissions\" dialog", () => grantSharingUserEdit(page, el("\"Edit Permissions\" dialog")));
    await session.step(44, "Then the sharing user's row of \"Edit Permissions\" dialog should show the user's picture", () => sharingUserPicture(page, el("\"Edit Permissions\" dialog")));
    await session.step(45, "When user clicks on CANCEL button in \"Edit Permissions\" dialog", () => clickOn(page, el("CANCEL button in \"Edit Permissions\" dialog")));
    await session.step(46, "Then \"Edit Permissions\" dialog should be absent", () => shouldBe(page, el("\"Edit Permissions\" dialog"), "absent"));
  });
});
