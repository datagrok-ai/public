/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/landing/pep-triage-campaigns.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {campaignsOnServer, cellIs, settingsDefaults, sortedBy, tableColumns} from '../../bindings/landing.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openApp, packageInstalled, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PepTriage landing: Hit Triage's campaigns table, and the table settings it shares with Hit Triage", () => {
  const session = feature(test, "features/landing/pep-triage-campaigns.feature", import.meta.url);
  test("PepTriage landing: Hit Triage's campaigns table, and the table settings it shares with Hit Triage", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the \"Bio\" package is installed", () => packageInstalled(page, "Bio"));
    await session.step(14, "And the campaigns tables start from their defaults, put back at feature end", () => settingsDefaults(page));
    await session.step(15, "And the \"PepTriage\" campaigns are on the server:", () => campaignsOnServer(page, "PepTriage", [["code","name","status","created","molecules","author","can edit"],[session.text("BDDPTL{time}-1"),session.text("BDD PT Alpha {time}"),session.text("BDD Shared {time}"),"2099/1/13","3","me","everyone"],[session.text("BDDPTL{time}-2"),session.text("BDD PT Locked {time}"),"In Progress","2099/1/12","5","sharing user","author"]]), [["code","name","status","created","molecules","author","can edit"],[session.text("BDDPTL{time}-1"),session.text("BDD PT Alpha {time}"),session.text("BDD Shared {time}"),"2099/1/13","3","me","everyone"],[session.text("BDDPTL{time}-2"),session.text("BDD PT Locked {time}"),"In Progress","2099/1/12","5","sharing user","author"]]);
    await session.step(19, "And the \"Hit Triage\" campaigns are on the server:", () => campaignsOnServer(page, "Hit Triage", [["code","name","status","created","molecules"],[session.text("BDDPTH{time}-1"),session.text("BDD PT Cross {time}"),session.text("BDD Shared {time}"),"2099/1/13","2"]]), [["code","name","status","created","molecules"],[session.text("BDDPTH{time}-1"),session.text("BDD PT Cross {time}"),session.text("BDD Shared {time}"),"2099/1/13","2"]]);
    await session.step(22, "And the campaigns table of \"PepTriage\" is sorted by \"Created\" descending", () => sortedBy(page, "PepTriage", "Created", "descending"));
    await session.step(23, "And user opens the \"PepTriage\" app", () => openApp(page, "PepTriage"));
    await run.scenario("The PepTriage landing lists its campaigns with Hit Triage's columns and row actions", async () => {
      await session.step(26, "Then the \"PepTriage\" view should be current", () => viewIsCurrent(page, "PepTriage"));
      await session.step(27, "And the campaigns table should have the columns \"Name, Created, Author, Total, Selected, Status\"", () => tableColumns(page, "Name, Created, Author, Total, Selected, Status"));
      await session.step(28, "And the \"Total\" column of the \"BDD PT Locked {time}\" campaign row should be \"5\"", () => cellIs(page, "Total", session.text("BDD PT Locked {time}"), "5"));
      await session.step(29, "And permissions icon of \"BDD PT Alpha {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("permissions icon of \"BDD PT Alpha {time}\" campaign row")), "visible"));
      await session.step(30, "And permissions icon of \"BDD PT Locked {time}\" campaign row should be hidden", () => shouldBe(page, el(session.text("permissions icon of \"BDD PT Locked {time}\" campaign row")), "hidden"));
    });
    await run.scenario("A grouping chosen in PepTriage groups the Hit Triage table as well", async () => {
      await session.step(33, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(34, "And user picks \"Status\" from the open menu", () => pickFromOpenMenu(page, "Status"));
      await session.step(35, "Then \"BDD Shared {time}\" campaign group should be visible", () => shouldBe(page, el(session.text("\"BDD Shared {time}\" campaign group")), "visible"));
      await session.step(36, "When user opens the \"Hit Triage\" app", () => openApp(page, "Hit Triage"));
      await session.step(37, "Then \"BDD Shared {time}\" campaign group should be visible", () => shouldBe(page, el(session.text("\"BDD Shared {time}\" campaign group")), "visible"));
      await session.step(38, "When user clicks on campaign grouping icon", () => clickOn(page, el("campaign grouping icon")));
      await session.step(39, "And user picks \"None\" from the open menu", () => pickFromOpenMenu(page, "None"));
      await session.step(40, "Then \"BDD PT Cross {time}\" campaign row should be visible", () => shouldBe(page, el(session.text("\"BDD PT Cross {time}\" campaign row")), "visible"));
      await session.step(41, "And \"BDD Shared {time}\" campaign group should be absent", () => shouldBe(page, el(session.text("\"BDD Shared {time}\" campaign group")), "absent"));
    });
    run.finish();
  });
});
