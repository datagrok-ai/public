# Login — Playwright Test Reference

**File:** `playwright-tests/e2e/login/login.test.ts`
**Target:** `https://dev.datagrok.ai` (configured in `playwright-tests/.env`)
**Total tests:** 27 (27 active, 0 skipped)

All tests run in an **unauthenticated** browser context — the global storageState from
`playwright.config.ts` is overridden with `{ cookies: [], origins: [] }` so each test
sees a fresh login page.

---

## Section 1 — Positive (5 tests)

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 1 | valid credentials show the ribbon | Fill `.env` credentials → click Login → `.d4-ribbon` visible within 30 s |
| 2 | pressing Enter in the password field submits the form | Fill credentials → press Enter in password field → `.d4-ribbon` visible |
| 3 | login page has Datagrok in its title | `page.title()` matches `/datagrok/i` |
| 4 | login form fields are empty on first load | Both input values are `""` on fresh navigation |
| 5 | login button is visible and enabled on page load | Login button is both visible and enabled before any interaction |

---

## Section 2 — Negative: empty and blank inputs (6 tests)

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 6 | clicking login with both fields empty does not authenticate | Click Login with no input → login form still visible, ribbon absent |
| 7 | valid login and empty password does not authenticate | Fill only username → click Login → not authenticated |
| 8 | empty login and valid password does not authenticate | Fill only password → click Login → not authenticated |
| 9 | whitespace-only login does not authenticate | Username = `"   "` + valid password → not authenticated |
| 10 | whitespace-only password does not authenticate | Valid username + password = `"   "` → not authenticated |
| 11 | both fields filled with whitespace only does not authenticate | Both `"   "` → not authenticated |

---

## Section 3 — Negative: wrong credentials (6 tests)

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 12 | wrong username and wrong password shows Login failed | Bogus credentials → error label contains "Login failed" |
| 13 | valid username with wrong password shows Login failed | `.env` username + wrong password → "Login failed" |
| 14 | wrong username with valid password shows Login failed | Bogus username + `.env` password → "Login failed" |
| 15 | correct username with password in wrong case shows Login failed | `.env` password uppercased → "Login failed" (password is case-sensitive) |
| 16 | numeric-only username shows Login failed | `"1234567890"` / `"1234567890"` → not authenticated |
| 17 | email-format username with wrong domain shows Login failed | `notreal@notexist.xyz` → not authenticated |

---

## Section 4 — Negative: boundary inputs (3 tests)

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 18 | very long username (500 chars) shows error without crash | 500-char username → no JS error, no redirect, login form still visible |
| 19 | very long password (500 chars) shows error without crash | 500-char password → no JS error, no redirect, login form still visible |
| 20 | single character username shows Login failed | `"x"` / `"y"` → not authenticated |

---

## Section 5 — Negative: special inputs (8 tests)

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 21 | SQL injection in login field shows error without crash | `' OR '1'='1` in both fields → not authenticated, no crash |
| 22 | SQL injection UNION SELECT in login field shows error without crash | `admin'--` in username → not authenticated, no crash |
| 23 | XSS attempt in login field does not execute script and shows error | `<script>alert(1)</script>` → no dialog fired, not authenticated |
| 24 | HTML injection in login field shows error without markup rendering | `<b>admin</b>` → not authenticated, no HTML rendered as markup |
| 25 | unicode characters in login field show Login failed | Cyrillic username/password → not authenticated |
| 26 | emoji in login field shows error without crash | `"😀user"` → no crash, not authenticated |
| 27 | null-byte-like string in login field shows error without crash | `"user\u0000admin"` → no crash, not authenticated |
| 28 | password with special characters shows Login failed for wrong user | `!@#$%^&*()` password with bogus user → not authenticated |

---

## Section 6 — Negative: repeated failures (2 tests)

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 29 | three consecutive failed logins keep the form accessible | 3 wrong logins in a row → login button still visible and enabled (no lockout / form destruction) |
| 30 | failed login followed by valid credentials succeeds | Wrong login first → error shown → re-fill with valid credentials → ribbon visible |

---

## Key CSS selectors

| Selector | Element |
|----------|---------|
| `#signup-login-fields input[placeholder="Login or Email"]` | Username / email input field |
| `#signup-login-fields input[placeholder="Password"]` | Password input field |
| `#signup-login-fields .signup-buttons button` | Login submit button |
| `#signup-login-fields label` | Error label (shows "Login failed" on auth failure) |
| `.d4-ribbon` | Main application ribbon — visible only after successful login |

---

## Coverage gaps (not automated)

| Scenario | Reason |
|----------|--------|
| "Forgot password" / password reset flow | Requires email access; out of scope for login-only tests |
| SSO / OAuth login | Requires third-party identity provider setup |
| Account lockout after N failed attempts | Lockout policy not confirmed; would require cleanup |
| "Remember me" / persistent session | No such checkbox observed in current UI |
| Mobile / responsive layout | Not covered by current test project (Desktop Chrome only) |
| Logout and redirect back to login | Covered implicitly by auth state; dedicated test pending |
| Admin-created account first login | Requires admin API; covered in registration-test.md |
