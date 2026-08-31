# npm trusted publishing (OIDC)

Everything this repo publishes to npm — `datagrok-api`, `datagrok-tools`, the
`@datagrok-libraries/*` libraries, the `@datagrok-misc/*` helpers and the
`@datagrok/*` plugins — authenticates to the registry with a short-lived
OpenID Connect token minted by GitHub Actions for that specific workflow run.

**There is no npm token in CI.** The workflows carry no `NODE_AUTH_TOKEN`, and
there is no `NPM_TOKEN` secret to leak. npm trusts *a workflow file in a
repository*, so a run of any other workflow — or of these workflows from a fork
— cannot publish.

Trusted publishing also turns on **provenance**: npm attests each tarball
against the run that produced it and shows a verified "Built and signed on
GitHub Actions" badge on npmjs.com. The npm CLI enables this automatically when
publishing over OIDC from a public repo; the workflows do not opt into it.

## Who publishes what

npm allows exactly **one** trusted publisher per package, so every package must
map to exactly one workflow. It already does:

| Source                | npm name(s)                   | Workflow file    |
|-----------------------|-------------------------------|------------------|
| `js-api/`             | `datagrok-api`                | `js-api.yml`     |
| `tools/`              | `datagrok-tools`              | `tools.yml`      |
| `libraries/<lib>/`    | `@datagrok-libraries/<lib>`   | `libraries.yaml` |
| `misc/<project>/`     | `@datagrok-misc/<project>`    | `misc.yaml`      |
| `packages/<Package>/` | `@datagrok/<package>`         | `packages.yaml`  |

The `libraries` and `packages` jobs in `js-api.yml` only bump the `datagrok-api`
dependency and commit — they never publish, so they need no configuration.

## What the workflows do

Each publishing job carries:

```yaml
    permissions:
      contents: write     # the job pushes the refreshed package-lock.json
      id-token: write     # lets the runner mint the OIDC token
```

and, immediately before `npm publish`, swaps the toolchain:

```yaml
      - name: Setup Node for trusted publishing
        uses: actions/setup-node@v4
        with:
          node-version: '24.x'
          registry-url: 'https://registry.npmjs.org'
      - run: npm install -g npm@latest
```

Trusted publishing needs **npm >= 11.5.1 on Node >= 22.14**, which is newer than
the Node the packages build against. Swapping only for the publish step keeps
every build on its pinned Node.

The publish steps themselves are plain `npm publish` with no `env:` block. The
npm CLI mints the OIDC token, exchanges it for a registry token scoped to that
one package, and publishes. A package that is not enrolled fails the publish
with a `401` — loudly, rather than silently falling back to anything.

## Enrolling a new package

npm cannot register a trusted publisher for a package that does not exist yet —
there is no equivalent of PyPI's pending publishers. So **the first version of a
brand-new package is published by hand**, then enrolled, and CI takes over from
the second version onwards.

Do this *before* merging to `master`. Every workflow skips publishing a version
that is already on the registry, so the merge that introduces the package is a
no-op for the publish step instead of a red build.

### 1. Point `repository` at this repo

npm's provenance check rejects a publish whose `repository.url` does not resolve
to the repository the build ran in (the comparison is case-sensitive). Every
`package.json` here needs:

```json
  "repository": {
    "type": "git",
    "url": "https://github.com/datagrok-ai/public.git"
  },
```

Optionally add `"directory": "packages/<Package>"` so npmjs.com links to the
right subfolder.

### 2. Publish the first version by hand

Set the version to `1.0.0` — `packages.yaml` refuses to publish anything below
`1.0.0`, so a package left at `0.x` never publishes from CI.

```bash
cd packages/<Package>
npm install
npm run build
npm login                # account with write access to the scope, 2FA on
npm publish --access public
```

This one publish has no provenance attached; every later version does.

### 3. Register the trusted publisher

```bash
npm install -g npm@latest        # npm trust needs >= 11.15.0
npm trust github @datagrok/<package> \
  --repo datagrok-ai/public \
  --file packages.yaml \
  --allow-publish --yes
```

`--file` is the **basename** of the workflow that publishes it — see the table
above. Use `libraries.yaml`, `misc.yaml`, `js-api.yml` or `tools.yml` for the
other source trees.

Or run the script, which resolves the workflow for you and skips anything
already enrolled:

```bash
./.github/scripts/npm-trust-enroll.sh @datagrok/<package>
```

### 4. Verify

```bash
npm trust list @datagrok/<package>
```

Then merge. On the next version bump, the job log shows `Publishing to
https://registry.npmjs.org with tag latest and public access` followed by a
Sigstore transparency-log URL, and npmjs.com shows the provenance badge.

## Bulk enrollment

**Run this from a real terminal on a machine with a browser.** Every `npm trust`
call is gated on a one-time password: npm prints a `https://www.npmjs.com/auth/cli/…`
URL and waits for you to approve it. It fails with `EOTP` in any non-interactive
shell, and a token that bypasses 2FA will not work either.

There are two equivalent scripts — note that the flag syntax differs:

```powershell
npm install -g npm@latest        # npm trust needs >= 11.15.0
npm login
.\.github\scripts\npm-trust-enroll.ps1 -DryRun    # review the plan
.\.github\scripts\npm-trust-enroll.ps1 | Tee-Object enroll.log
```

```bash
npm install -g npm@latest
npm login
./.github/scripts/npm-trust-enroll.sh --dry-run   # review the plan
./.github/scripts/npm-trust-enroll.sh | tee enroll.log
```

Both walk `js-api`, `tools`, `libraries/*`, `misc/*` and `packages/*`, skip
packages that are private, not Datagrok-owned, or not yet on the registry, and
sleep 2s between writes to stay under the rate limit. The dry run never contacts
the registry, so it prints the plan but proves nothing about auth.

Exit codes: `0` all good · `1` some packages failed · `2` `-Only`/filter matched
nothing · `3` stopped because the first package failed with nothing yet
succeeded (systemic — almost always the OTP gate).

On the first OTP prompt npmjs.com offers **"skip two-factor authentication for
the next 5 minutes"** — enable it and the rest of the run goes through
unattended (roughly 80 packages fit in the window). If it expires mid-run, just
re-run: packages already configured come back as `SKIP … already has a trusted
publisher`.

Nothing publishes until enrollment is done, so run this before the next release
train. Anything still missing shows up as a failed publish step, not a silent
one.

## Optional hardening: a GitHub environment

`npm trust github` also accepts `--env <name>`. Add `environment: <name>` to the
publishing jobs and the same name to the trusted publisher, and npm will only
accept publishes from runs that passed that environment's protection rules
(branch restrictions, required reviewers, wait timers). Not enabled today: the
jobs already gate publishing on `github.ref == 'refs/heads/master'`, and pull
requests from forks never receive `id-token: write`.

Worth doing at the npm org level as well: set **"Require two-factor
authentication and disallow tokens"** so no one can mint a publishing token that
bypasses this.

## Troubleshooting

| Symptom | Cause |
|---|---|
| `401 Unauthorized` on publish | The OIDC exchange did not produce a token — usually the package has no trusted publisher, or its publisher points at a different workflow file. Check `npm trust list <pkg>`. |
| `ENEEDAUTH` / `This command requires you to be logged in` | Something defined `NODE_AUTH_TOKEN` as an empty string. An *undefined* `NODE_AUTH_TOKEN` is fine; an empty one makes npm think auth is configured and skips the OIDC token. |
| `422 Unprocessable Entity` mentioning provenance or the repository | `repository.url` in `package.json` does not match `datagrok-ai/public`. |
| `Provenance generation in GitHub Actions requires "write" access to the "id-token" permission` | The `permissions:` block on the job is missing `id-token: write`. |
| `npm trust` says a configuration already exists | One publisher per package. `npm trust list <pkg>`, then `npm trust revoke <pkg> --id <id>`, then re-add. |
| `npm trust` fails with a 2FA error | It needs account-level 2FA; granular access tokens with "bypass 2FA" are rejected. |

Add `--loglevel verbose` to `npm publish` to see the `oidc` lines explaining
exactly which step failed.

## References

* [npm trusted publishers](https://docs.npmjs.com/trusted-publishers)
* [`npm trust`](https://docs.npmjs.com/cli/v11/commands/npm-trust)
* [Generating provenance statements](https://docs.npmjs.com/generating-provenance-statements)
