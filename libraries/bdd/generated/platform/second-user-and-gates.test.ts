/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/platform/second-user-and-gates.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {packageInstalled, runningAccountSignedIn, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, standRunsService, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The second account and the capability gates", () => {
  const session = feature(test, "features/platform/second-user-and-gates.feature", import.meta.url);
  test("The sharing user signs in on the same page and the running account comes back", {tag: ["@platform"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "When user signs in as the sharing user", () => signInAsSharingUser(page));
    await session.step(12, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
    await session.step(13, "And no errors should have been logged", () => noErrors(page));
    await session.step(14, "When user signs in as themselves again", () => signInAsSelf(page));
    await session.step(15, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    await session.step(16, "And no errors should have been logged", () => noErrors(page));
  });
  test("A service the stand does not run skips the rest of the test", {tag: ["@platform"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "Given the stand runs the \"No Such Service\" service", () => standRunsService(page, "No Such Service"));
    await session.step(20, "Then the \"No such view\" view should be current", () => viewIsCurrent(page, "No such view"));
  });
  test("A package the stand does not have skips the rest of the test", {tag: ["@platform"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "Given the \"NoSuchPackageAnywhere\" package is installed", () => packageInstalled(page, "NoSuchPackageAnywhere"));
    await session.step(24, "Then the \"No such view\" view should be current", () => viewIsCurrent(page, "No such view"));
  });
});
