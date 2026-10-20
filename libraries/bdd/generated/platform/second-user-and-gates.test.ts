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
import {packageInstalled, runningAccountSignedIn, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The second account and the capability gates", () => {
  const session = feature(test, "features/platform/second-user-and-gates.feature", import.meta.url);
  test("The sharing user signs in on the same page and the running account comes back", {tag: ["@platform"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "When user signs in as the sharing user", () => signInAsSharingUser(page));
    await session.step(15, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
    await session.step(16, "And no errors should have been logged", () => noErrors(page));
    await session.step(17, "When user signs in as themselves again", () => signInAsSelf(page));
    await session.step(18, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    await session.step(19, "And no errors should have been logged", () => noErrors(page));
  });
  test("A package the stand does not have skips the rest of the test", {tag: ["@platform"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "Given the \"NoSuchPackageAnywhere\" package is installed", () => packageInstalled(page, "NoSuchPackageAnywhere"));
    await session.step(23, "Then the \"No such view\" view should be current", () => viewIsCurrent(page, "No such view"));
  });
});
