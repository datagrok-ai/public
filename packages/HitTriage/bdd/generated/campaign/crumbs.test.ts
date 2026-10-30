/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/campaign/crumbs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/landing.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {backByCrumb, hitDesignCampaign, openByLink} from '../../bindings/campaign.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {urlShouldContain, urlShouldNotContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The crumbs of a Hit Design campaign view go back to the landing", () => {
  const session = feature(test, "features/campaign/crumbs.feature", import.meta.url);
  test("The first crumb goes back to the landing and closes the campaign", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the Hit Design campaign \"BDD Crumbs\" with the stages \"Idea\" and a table of:", () => hitDesignCampaign(page, "BDD Crumbs", "Idea", [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"]]), [["Molecule","Stage","V-iD","MW"],["CCO","Idea","","46.07"]]);
    await session.step(15, "When user opens the campaign \"BDD Crumbs\" by its link", () => openByLink(page, "BDD Crumbs"));
    await session.step(18, "Then the page address should contain \"campaignId=BDDBDDCRUMBS-1\"", () => urlShouldContain(page, "campaignId=BDDBDDCRUMBS-1"));
    await session.step(19, "When user goes back by the \"Hit Design\" crumb of the campaign \"BDD Crumbs\"", () => backByCrumb(page, "Hit Design", "BDD Crumbs"));
    await session.step(20, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
    await session.step(21, "And the open table views should be exactly \"\"", () => openTableViewsExactly(page, ""));
    await session.step(22, "And the page address should not contain \"campaignId\"", () => urlShouldNotContain(page, "campaignId"));
  });
});
