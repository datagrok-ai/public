---
title: "Back up and restore"
sidebar_position: 4
description: What to back up on a Datagrok instance, how to keep the backups consistent, and in which order to restore them.
keywords:
  - backup
  - restore
  - disaster recovery
  - postgres backup
  - server keys
  - file storage
---

A Datagrok instance keeps its state in three places. Back up all three, as
of the same point in time, to be able to restore the instance.

| What                   | Where it lives                                                                                          | Contains                                                                                  |
|------------------------|---------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------|
| **Datagrok database**  | The Postgres database the server connects to (Amazon RDS, Cloud SQL, or the bundled Postgres container) | Users, groups, permissions, connections, queries, dashboards, settings, audit log, and the encrypted credentials |
| **File storage**       | The S3 or Google Cloud Storage bucket, or the local data volume when no bucket is configured             | Uploaded files, table data of saved projects, package files, and cached results           |
| **Server keys**        | PEM files, such as `auth.pem` and `datagrok.pem`, in the `settings` folder of the bucket, or on the configuration volume (`datagrok-cfg` in the Helm chart) when no bucket is configured. Alternatively, AWS Secrets Manager or GCP Secret Manager if you [moved the keys](../../govern/access-control/server-keys.md#move-a-key-between-backends) there | The keys that encrypt stored credentials and sign sign-in tokens |

:::caution

The credentials in the database are encrypted with the server keys. A restored
database without the matching keys still works, but every stored password and
token is lost and must be entered again. Keep the keys backed up, and protected
at least as well as the database.

:::

Container images and their configuration are not state. Record which version
you run (see [Images and versions](../images.md)) so you can deploy the same
images during a restore.

## Consistency

Projects in the database point at table data in the file storage. A database
backup newer than the file storage backup contains projects whose data is
missing. To keep them in step:

* Take the database backup and the file storage backup on the same schedule,
  database first.
* For a guaranteed consistent pair, for example before an
  [upgrade](../upgrade.md), stop the Datagrok server containers, back up both,
  and start them again.
* Back up the whole database, not individual schemas. Plugins keep their own
  tables in it.

## By platform

| Platform                        | Database                                            | File storage and keys                              |
|---------------------------------|-----------------------------------------------------|----------------------------------------------------|
| AWS                             | RDS automated backups or snapshots                  | S3 versioning and AWS Backup ([setup](configure-s3-backup.md)) |
| Google Cloud                    | Cloud SQL automated backups                         | GCS object versioning                              |
| Kubernetes with the Helm chart  | `pg_dump`, or snapshots of the Postgres volume       | Snapshots of the `datagrok-data` and `datagrok-cfg` volumes |
| Docker Compose                  | `pg_dump` from the Postgres container               | A copy of the Datagrok data volume                 |

## Restore

1. Deploy the same Datagrok version the backup was taken from, with the server
   stopped.
1. Restore the database.
1. Restore the file storage and the server keys. If the keys live in a secrets
   manager, make sure the same secrets are available.
1. Start the server and sign in as an administrator.
1. Check `/settings/initial/health`, open a few dashboards, and test a data
   connection that uses stored credentials.

See also:

* [Upgrade and roll back](../upgrade.md)
* [Server keys](../../govern/access-control/server-keys.md)
* [Deployment](../deploy.md)
