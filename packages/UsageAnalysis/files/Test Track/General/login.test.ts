import { test, expect, Page } from '@playwright/test';
import * as dotenv from 'dotenv';

dotenv.config();

// Login tests run against a fresh (unauthenticated) browser context.
// This overrides the global storageState set in playwright.config.ts.
test.use({ storageState: { cookies: [], origins: [] } });

const BASE_URL = process.env.DATAGROK_URL!;
const VALID_LOGIN = process.env.DATAGROK_LOGIN!;
const VALID_PASSWORD = process.env.DATAGROK_PASSWORD!;

const SEL = {
  loginField: '#signup-login-fields input[placeholder="Login or Email"]',
  passwordField: '#signup-login-fields input[placeholder="Password"]',
  loginButton: '#signup-login-fields .signup-buttons button',
  // .signup-status is always in DOM (display:block) — empty text on load, "Login failed" after error
  errorLabel: '#signup-login-fields .signup-status',
  ribbon: '.d4-ribbon',
};

// ---------------------------------------------------------------------------
// Page-level helpers
// ---------------------------------------------------------------------------

async function gotoLogin(page: Page): Promise<void> {
  await page.goto(BASE_URL);
  await page.waitForSelector(SEL.loginField, { timeout: 15_000 });
}

async function fillAndSubmit(page: Page, username: string, password: string): Promise<void> {
  await page.locator(SEL.loginField).fill(username);
  await page.locator(SEL.passwordField).fill(password);
  await page.locator(SEL.loginButton).click();
}

/** Wait for ".signup-status" to show "Login failed" text (server round-trip). */
async function expectLoginFailed(page: Page): Promise<void> {
  await expect(page.locator(SEL.errorLabel)).toHaveText(/login failed/i, { timeout: 8_000 });
}

/** Confirm the user is NOT authenticated: login form still present, ribbon absent. */
async function expectNotLoggedIn(page: Page): Promise<void> {
  await expect(page.locator(SEL.loginField)).toBeVisible({ timeout: 5_000 });
  await expect(page.locator(SEL.ribbon)).not.toBeVisible();
}

// ---------------------------------------------------------------------------
// Positive scenarios
// ---------------------------------------------------------------------------

test.describe('Login — Positive', () => {
  test('valid credentials show the ribbon', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, VALID_LOGIN, VALID_PASSWORD);
    await expect(page.locator(SEL.ribbon).first()).toBeVisible({ timeout: 30_000 });
  });

  test('pressing Enter in the password field submits the form', async ({ page }) => {
    await gotoLogin(page);
    await page.locator(SEL.loginField).fill(VALID_LOGIN);
    await page.locator(SEL.passwordField).fill(VALID_PASSWORD);
    await page.locator(SEL.passwordField).press('Enter');
    await expect(page.locator(SEL.ribbon).first()).toBeVisible({ timeout: 30_000 });
  });

  test('login page has Datagrok in its title', async ({ page }) => {
    await gotoLogin(page);
    await expect(page).toHaveTitle(/datagrok/i);
  });

  test('login form fields are empty on first load', async ({ page }) => {
    await gotoLogin(page);
    await expect(page.locator(SEL.loginField)).toHaveValue('');
    await expect(page.locator(SEL.passwordField)).toHaveValue('');
  });

  test('login button is visible and enabled on page load', async ({ page }) => {
    await gotoLogin(page);
    await expect(page.locator(SEL.loginButton)).toBeVisible();
    await expect(page.locator(SEL.loginButton)).toBeEnabled();
  });

  test('error label is empty on fresh page load', async ({ page }) => {
    await gotoLogin(page);
    await expect(page.locator(SEL.errorLabel)).toHaveText('');
  });
});

// ---------------------------------------------------------------------------
// Negative: empty / blank inputs
// ---------------------------------------------------------------------------

test.describe('Login — Negative: empty and blank inputs', () => {
  test('clicking login with both fields empty shows Login failed', async ({ page }) => {
    await gotoLogin(page);
    await page.locator(SEL.loginButton).click();
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });

  test('valid login and empty password shows Login failed', async ({ page }) => {
    await gotoLogin(page);
    await page.locator(SEL.loginField).fill(VALID_LOGIN);
    await page.locator(SEL.loginButton).click();
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });

  test('empty login and valid password shows Login failed', async ({ page }) => {
    await gotoLogin(page);
    await page.locator(SEL.passwordField).fill(VALID_PASSWORD);
    await page.locator(SEL.loginButton).click();
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });

  test('whitespace-only login does not authenticate', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, '   ', VALID_PASSWORD);
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });

  test('whitespace-only password does not authenticate', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, VALID_LOGIN, '   ');
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });

  test('both fields filled with whitespace only does not authenticate', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, '   ', '   ');
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });
});

// ---------------------------------------------------------------------------
// Negative: wrong credentials
// ---------------------------------------------------------------------------

test.describe('Login — Negative: wrong credentials', () => {
  test('wrong username and wrong password shows Login failed', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, 'no_such_user_xyz', 'wrongpass123');
    await expectLoginFailed(page);
  });

  test('valid username with wrong password shows Login failed', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, VALID_LOGIN, 'wrongpass_123!');
    await expectLoginFailed(page);
  });

  test('wrong username with valid password shows Login failed', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, 'no_such_user_xyz', VALID_PASSWORD);
    await expectLoginFailed(page);
  });

  test('correct username with password in wrong case shows Login failed', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, VALID_LOGIN, VALID_PASSWORD.toUpperCase());
    await expectLoginFailed(page);
  });

  test('numeric-only username shows Login failed', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, '1234567890', '1234567890');
    await expectLoginFailed(page);
  });

  test('email-format username with wrong domain shows Login failed', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, 'notreal@notexist.xyz', 'somepassword');
    await expectLoginFailed(page);
  });
});

// ---------------------------------------------------------------------------
// Negative: boundary / extreme inputs
// ---------------------------------------------------------------------------

test.describe('Login — Negative: boundary inputs', () => {
  test('very long username (500 chars) shows error without crash', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, 'a'.repeat(500), 'password');
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });

  test('very long password (500 chars) shows error without crash', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, 'user', 'p'.repeat(500));
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });

  test('single character credentials show Login failed', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, 'x', 'y');
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });
});

// ---------------------------------------------------------------------------
// Negative: special / potentially dangerous inputs
// ---------------------------------------------------------------------------

test.describe('Login — Negative: special inputs', () => {
  test('SQL injection in login field shows error without crash', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, "' OR '1'='1", "' OR '1'='1");
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });

  test('SQL comment injection in login field shows error without crash', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, "admin'--", 'anypassword');
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });

  test('XSS attempt in login field does not execute script and shows error', async ({ page }) => {
    let alertFired = false;
    page.on('dialog', async (dialog) => {
      alertFired = true;
      await dialog.accept();
    });
    await gotoLogin(page);
    await fillAndSubmit(page, '<script>alert(1)</script>', 'password');
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
    expect(alertFired).toBe(false);
  });

  test('HTML tag in login field shows error without rendering markup', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, '<b>admin</b>', 'password');
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });

  test('unicode / cyrillic characters in login field show Login failed', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, 'пользователь', 'пароль');
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });

  test('emoji in login field shows error without crash', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, '😀user', 'password');
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });

  test('null-byte character in login field shows error without crash', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, 'user\u0000admin', 'password');
    // Server may respond with "Login failed" or "Operation caused an exception" — either is acceptable
    await expect(page.locator(SEL.errorLabel)).not.toHaveText('', { timeout: 8_000 });
    await expectNotLoggedIn(page);
  });

  test('special characters in password show Login failed', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, 'no_such_user', '!@#$%^&*()_+-=[]{}|;:,.<>?');
    await expectLoginFailed(page);
    await expectNotLoggedIn(page);
  });
});

// ---------------------------------------------------------------------------
// Negative: repeated / sequential failures
// ---------------------------------------------------------------------------

test.describe('Login — Negative: repeated failures', () => {
  test('three consecutive failed logins keep the form accessible', async ({ page }) => {
    await gotoLogin(page);
    for (let i = 0; i < 3; i++) {
      await page.locator(SEL.loginField).fill(`wrong_user_${i}`);
      await page.locator(SEL.passwordField).fill(`wrong_pass_${i}`);
      await page.locator(SEL.loginButton).click();
      await expect(page.locator(SEL.errorLabel)).toHaveText(/login failed/i, { timeout: 8_000 });
    }
    await expectNotLoggedIn(page);
    await expect(page.locator(SEL.loginButton)).toBeVisible();
    await expect(page.locator(SEL.loginButton)).toBeEnabled();
  });

  test('failed login followed by valid credentials succeeds', async ({ page }) => {
    await gotoLogin(page);
    await fillAndSubmit(page, 'no_such_user_xyz', 'wrongpass');
    await expect(page.locator(SEL.errorLabel)).toHaveText(/login failed/i, { timeout: 8_000 });
    await page.locator(SEL.loginField).fill(VALID_LOGIN);
    await page.locator(SEL.passwordField).fill(VALID_PASSWORD);
    await page.locator(SEL.loginButton).click();
    await expect(page.locator(SEL.ribbon).first()).toBeVisible({ timeout: 30_000 });
  });
});
