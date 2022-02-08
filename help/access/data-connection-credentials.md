<!-- TITLE: Data connection credentials-->
<!-- SUBTITLE: -->

# Credentials

Any used data source requires authentication to access the data. Credentials (logins, passwords, etc) can be specified
directly (Manual) or through integrating with secret manager services storing these credentials.

## Option "Manual"

This basic option provides keeping encrypted credentials for any data sources in the protected service inside the
Datagrok platform.

## Secrets Managers

Other options suggest using API calls to cloud secrets managers for authentication data sources. The primary purpose of
a secret management services is to allow you to decouple storage of secrets from the code or configuration that consumes
the secrets. This decoupling supports centralization, revocation and rotation passwords.

The Datagrok platform is integrated with storing secrets services:

* Amazon AWS Secrets Manager

The platform supports having several Secrets Managers Connection used by different groups. In these cases, specific
credentials for Data Source is stored in Cloud Manager that is referenced by selected Datagrok connection and specified
by Secret Name. The specified Secret Name should be defined in AWS Secrets Manager. The correct name of the Secret can
be evaluated during the creation of a connection by testing.

## How to setup connection using AWS Secret Manager

There is a Postgres database biodb01 located somewhere and available by credentials user01/password01

Secret with Name `biodb01` is created in AWS Secrets Manager using command:

```aws secretsmanager create-secret --name biodb01 --secret-string "{\"password\":\"password01\",\"login\":\"user01\"}"```

In Datagrok Platform we should perform the following actions:

1. Create connection aws-sm-bio to AWS Secret Manager (later we can use this connection for other secret names)
   ![Create connection to AWS Secret Manager](data-connection-secret-p01.png)

2. Create connection type Postgres to database biodb01 with name biodb01-connection. It is important to select
   aws-sm-bio created at the previous step and specify the correct secret name.  
   ![Add new Data connection](data-connection-secret-p02.png)

3. Click test to ensure that secret name is correct.
   ![Test connection](data-connection-secret-p03.png)

4. Click OK. Save setting.
