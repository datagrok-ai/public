---
title: "Secrets Managers"
format: mdx
sidebar_position: 2
---

# Secrets Managers

Datagrok supports secure connections to **AWS** and **GCP**, allowing you to:

- Access cloud services such as **BigQuery**, **S3**, and **Athena**.
- Configure log exports to **CloudWatch** or **Google Cloud Logging**.
- Securely provide credentials for other Datagrok connections by fetching them from the respective **Secrets Manager**.

> **Note:** You must be a member of the **Administrators** or **Developers** group to create these connections.

---

## AWS Connection

Datagrok supports two authentication methods for AWS:

1. **IAM credentials** – Manually provide an **Access Key ID** and **Secret Access Key**.
2. **Task role credentials** – When Datagrok is running on **EC2** or **ECS/Fargate**, it can fetch temporary credentials automatically from the IAM role assigned to the environment.
    - **EC2:** Credentials are retrieved from the Instance Metadata Service (IMDS).
    - **ECS/Fargate:** Credentials are fetched from the task role endpoint.

### Creating an AWS Connection
1. Navigate to **Platform → Credentials → AWS**.
2. Right-click the provider and select **Add connection...**.
3. Select the **Authentication Method**.
4. Fill in the required fields and **save** the connection.

### Using an AWS Connection
The IAM or task role used by Datagrok must have the necessary permissions for the resources you intend to use, such as:

- **Athena:** `athena:*`, `s3:*`, and `glue:*` (for the Data Catalog).
- **Logging:** `logs:CreateLogGroup`, `logs:CreateLogStream`, `logs:PutLogEvents`.

#### As credentials provider for Athena
1. Navigate to **Databases → Athena**.
2. Right-click the provider and select **Add connection...**.
3. Set **Credentials** to the AWS connection you created, then complete the remaining parameters.
4. (Optional) Specify a **Secret name**. Datagrok will fetch the credentials object from [AWS Secrets Manager](https://docs.aws.amazon.com/secretsmanager/latest/userguide/introduction.html) using this **AWS** connection.

#### For logs export to CloudWatch
1. Navigate to **Settings → Log → Log Export → Add New Export Block**.
2. Select **Amazon CloudWatch** and choose the AWS connection in the **Connection** field.

---

## GCP Connection

Datagrok supports two authentication methods for GCP:

1. **Service account key** – Upload a JSON key file for a GCP service account.
2. **Service account impersonation** – Provide a **Service Account Email**. If Datagrok runs under a service account with the **Service Account Token Creator** role, it can impersonate the specified account.

### Creating a GCP Connection
1. Navigate to **Databases → GCP**.
2. Right-click the provider and select **Add connection...**.
3. Select the **Authentication Method**.
4. Fill in the required fields and **save** the connection.

### Using a GCP Connection
The service account (or impersonated account) must have the necessary permissions:

- **BigQuery:** `roles/bigquery.dataViewer` (dataset access).
- **Logging:** `roles/logging.logWriter` (permission to write logs).

#### As credentials provider for BigQuery
1. Navigate to **Databases → BigQuery**.
2. Right-click the provider and select **Add connection...**.
3. Set **Credentials** to the GCP connection you created, then complete the remaining parameters.
4. (Optional) Specify a **Secret name**. Datagrok will fetch the credentials object from [Google Secret Manager](https://cloud.google.com/security/products/secret-manager) using this **GCP** connection.

#### For logs export to Google Cloud Logging
1. Navigate to **Settings → Log → Log Export → Add New Export Block**.
2. Select **Google Cloud Logging** and choose the GCP connection in the **Connection** field.

---

See also:

* [Data connection](../../access/access.md#data-connection)
* [Adding connection](../../access/databases/databases.md#adding-connection)
