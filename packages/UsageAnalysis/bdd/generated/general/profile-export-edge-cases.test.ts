/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/profile-export-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, clickOn, downloadContains, downloadThrough, enterInto, fileDownloaded, shouldBe, uncheck, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {dialogCloses, openTableOf, packageInstalled, runningAccountSignedIn, signInAs, signInAsSelf, userNameOnServer, userOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Profile and export edge cases: an empty and a long name, an empty molecule in a CSV", () => {
  const session = feature(test, "features/general/profile-export-edge-cases.feature", import.meta.url);
  test("An emptied first name keeps OK disabled", {tag: ["@serial", "@edge-cases", "@known-failure:27", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(20, "Given a user \"bddprofile\" is on the server", () => userOnServer(page, "bddprofile"));
      await session.step(21, "When user signs in as \"bddprofile\"", () => signInAs(page, "bddprofile"));
      await session.step(22, "And user opens the address \"/u\"", () => openAddress(page, "/u"));
      await session.step(23, "Then \"Logout\" link should be visible", () => shouldBe(page, el("\"Logout\" link"), "visible"));
      await session.step(24, "When user clicks on second \"Edit\" icon", () => clickOn(page, el("second \"Edit\" icon")));
      await session.step(25, "Then \"Enter new name\" dialog should be visible", () => shouldBe(page, el("\"Enter new name\" dialog"), "visible"));
      await session.step(26, "When user clears \"First name\" input in \"Enter new name\" dialog", () => clearField(page, el("\"First name\" input in \"Enter new name\" dialog")));
      await session.step(27, "Then OK button in \"Enter new name\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Enter new name\" dialog"), "disabled"));
      await session.step(28, "When user clicks on CANCEL button in \"Enter new name\" dialog", () => clickOn(page, el("CANCEL button in \"Enter new name\" dialog")));
      await session.step(29, "Then the \"Enter new name\" dialog should close", () => dialogCloses(page, "Enter new name"));
      await session.step(30, "And the user \"bddprofile\" should have the name \"bddprofile\" on the server", () => userNameOnServer(page, "bddprofile", "bddprofile"));
      await session.step(31, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(32, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    }, 27);
  });
  test("A 200-character first name is saved whole", {tag: ["@serial", "@edge-cases", "@known-failure:46", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(38, "Given a user \"bddprofile\" is on the server", () => userOnServer(page, "bddprofile"));
      await session.step(39, "When user signs in as \"bddprofile\"", () => signInAs(page, "bddprofile"));
      await session.step(40, "And user opens the address \"/u\"", () => openAddress(page, "/u"));
      await session.step(41, "Then \"Logout\" link should be visible", () => shouldBe(page, el("\"Logout\" link"), "visible"));
      await session.step(42, "When user clicks on second \"Edit\" icon", () => clickOn(page, el("second \"Edit\" icon")));
      await session.step(43, "Then \"Enter new name\" dialog should be visible", () => shouldBe(page, el("\"Enter new name\" dialog"), "visible"));
      await session.step(44, "When user enters \"LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL\" into \"First name\" input in \"Enter new name\" dialog", () => enterInto(page, "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL", el("\"First name\" input in \"Enter new name\" dialog")));
      await session.step(45, "And user clicks on OK button in \"Enter new name\" dialog", () => clickOn(page, el("OK button in \"Enter new name\" dialog")));
      await session.step(46, "Then the \"Enter new name\" dialog should close", () => dialogCloses(page, "Enter new name"));
      await session.step(47, "And the user \"bddprofile\" should have the name \"LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL\" on the server", () => userNameOnServer(page, "bddprofile", "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL"));
      await session.step(48, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(49, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(50, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    }, 46);
  });
  test("An empty molecule cell is exported as an empty field", {tag: ["@serial", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(53, "Given the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(54, "And user opens a table \"mols\" with:", () => openTableOf(page, "mols", [["smiles","id"],["CCO","1"],["","2"],["c1ccccc1","3"]]), [["smiles","id"],["CCO","1"],["","2"],["c1ccccc1","3"]]);
    await session.step(59, "Then \"smiles\" column should have semantic type \"Molecule\"", () => columnSemType(page, "smiles", "Molecule"));
    await session.step(60, "Given user watches downloads", () => watchDownloads(page));
    await session.step(61, "When user clicks on Export icon in toolbar", () => clickOn(page, el("Export icon in toolbar")));
    await session.step(62, "And user clicks on \"As CSV (options)...\" text in toolbar", () => clickOn(page, el("\"As CSV (options)...\" text in toolbar")));
    await session.step(63, "Then \"Save as CSV\" dialog should be visible", () => shouldBe(page, el("\"Save as CSV\" dialog"), "visible"));
    await session.step(64, "When user checks \"Molecules as Smiles\" checkbox in \"Save as CSV\" dialog", () => check(page, el("\"Molecules as Smiles\" checkbox in \"Save as CSV\" dialog")));
    await session.step(65, "And user unchecks \"Selected Columns Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Selected Columns Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(66, "And user unchecks \"Selected Rows Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Selected Rows Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(67, "And user unchecks \"Filtered Rows Only\" checkbox in \"Save as CSV\" dialog", () => uncheck(page, el("\"Filtered Rows Only\" checkbox in \"Save as CSV\" dialog")));
    await session.step(68, "And user downloads a file through OK button in \"Save as CSV\" dialog", () => downloadThrough(page, el("OK button in \"Save as CSV\" dialog")));
    await session.step(69, "Then a file \"mols.csv\" should have been downloaded", () => fileDownloaded(page, "mols.csv"));
    await session.step(70, "And the downloaded file \"mols.csv\" should contain text \",2\"", () => downloadContains(page, "mols.csv", ",2"));
    await session.step(71, "And the downloaded file \"mols.csv\" should contain text \"CCO,1\"", () => downloadContains(page, "mols.csv", "CCO,1"));
    await session.step(72, "And no errors should have been logged", () => noErrors(page));
  });
});
