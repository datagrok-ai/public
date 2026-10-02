---
title: "Upgrade and roll back"
sidebar_position: 10
description: How to plan, run, verify, and roll back a Datagrok upgrade on a self-hosted instance.
keywords:
  - upgrade datagrok
  - update datagrok version
  - rollback
  - downgrade
  - database migration
  - release planning
---

An upgrade replaces the Datagrok container images, and the configuration that
comes with them, and migrates the database on the first start. Database
migrations only go forward, so the way back is a restore from a backup. Plan
every production upgrade with that in mind.

## Before you upgrade

1. **Read the release notes.** Check [What's new](releases/releases.md) and the
   [release history](releases/release-history.md) for the target version. To
   check the JS API changes between two versions that can affect your own
   plugins, use the [compatibility tool](releases/compatibility/compatibility.mdx).
1. **Get the configuration from the same release.** A Compose file, Helm chart,
   or CloudFormation template written for one release can be invalid for
   another. See
   [Configuration is versioned with the image](images.md#release-binding).
1. **Go one minor version at a time.** For production, upgrade through
   consecutive minor versions, for example 1.26 → 1.27 → 1.28.
1. **Try it on a validation instance.** Restore a recent production backup to a
   separate instance, upgrade it, and run your checks (see
   [After the upgrade](#after-the-upgrade)).
1. **Back up production.** Take a consistent backup of the database, the file
   storage, and the server keys right before the upgrade. See
   [Back up and restore](complete-setup/backup.md).

## Upgrade

| Deployment              | How                                                                                                  |
|-------------------------|------------------------------------------------------------------------------------------------------|
| Helm chart              | `helm upgrade` with the chart tag `<version>-helm`. See [Upgrades](k8s/install-helm-chart.md#upgrades) |
| AWS CloudFormation (EKS) | Update the stack with the template for the new release. See [Update Datagrok components](aws/deploy-amazon-eks.mdx#update-datagrok-components) |
| Docker Compose          | Take the Compose file from the `release/<version>` branch, then `docker compose pull` and `docker compose up -d`. See [Images and versions](images.md) |

Upgrades don't replace the database or the file storage. Users, content,
credentials, and settings stay in place.

On the first start of the new version, the server applies the pending
database migrations before it accepts connections, and runs the remaining data
migrations right after startup. Expect the first start to take longer than
usual. Some migrations stop other running server instances first, so plan a
maintenance window rather than a rolling upgrade across versions.

Plugins installed with the `latest` version keep updating themselves from the
package registry, independently of platform upgrades. On production, consider
pinning plugin versions and updating them deliberately:
`grok s packages outdated` lists what's behind, and
`grok s packages install <name> --version <version>` pins a version (see
[Manage packages](../develop/server-management.md#manage-packages)).

## After the upgrade

1. Open `/settings/initial/health` on your instance and make sure every
   service is running, or
   run `grok s healthcheck` (see
   [Server health](../develop/server-management.md#server-health)).
1. Sign in through your identity provider, not only as `admin`.
1. Open the dashboards your organization depends on, and test their data
   connections.
1. Run your package tests (see
   [Test packages](../develop/how-to/tests/test-packages.md)).
1. Watch **Usage Analysis** > **Errors** for new errors in the following days.

## Roll back

Datagrok doesn't support running an older version against a database that a
newer version has migrated. To roll back:

1. Stop the Datagrok server containers.
1. Restore the database and the file storage from the backup you took before
   the upgrade.
1. Deploy the previous version with its own configuration (Compose file, chart
   tag, or template).
1. Start the server and repeat the checks from [After the upgrade](#after-the-upgrade).

Anything users created or changed after the upgrade is lost in a rollback.
Decide early: the longer production runs on the new version, the more a
rollback costs.

See also:

* [Deployment](deploy.md)
* [Images and versions](images.md)
* [Back up and restore](complete-setup/backup.md)
