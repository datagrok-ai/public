# Snowflake OAuth (JWT) Authentication

Datagrok supports per-user Snowflake access using OAuth (JWT) authentication. When enabled, each user's queries run under their own Snowflake identity, so row-level security policies and role-based access apply automatically.

## How It Works

1. User logs into Datagrok (via SSO, username/password, or any other method)
2. User runs a query against a Snowflake connection configured with OAuth (JWT)
3. Datagrok server mints a short-lived JWT containing the user's email, signed with a private key
4. Snowflake validates the JWT using the matching public key and maps the email to a Snowflake user
5. The query executes under that user's identity — row access policies, roles, and grants apply

## Prerequisites

- Datagrok instance with users that have email addresses configured
- Snowflake account with `ACCOUNTADMIN` access (for initial setup)
- Snowflake users with `LOGIN_NAME` matching Datagrok user emails

## Setup

### Step 1: Generate a Key Pair

Generate an RSA key pair on any machine:

```bash
openssl genrsa 2048 > datagrok_snowflake.pem
openssl rsa -in datagrok_snowflake.pem -pubout
```

Save the private key file (`datagrok_snowflake.pem`) — you will upload it to Datagrok.
Copy the public key output (the content between `-----BEGIN PUBLIC KEY-----` and `-----END PUBLIC KEY-----`) — you will paste it into Snowflake.

### Step 2: Configure Snowflake

Log into Snowflake as `ACCOUNTADMIN` and run the following, replacing the placeholders:

```sql
CREATE SECURITY INTEGRATION datagrok_jwt_oauth
  TYPE = EXTERNAL_OAUTH
  ENABLED = TRUE
  EXTERNAL_OAUTH_TYPE = CUSTOM
  EXTERNAL_OAUTH_ISSUER = '<your-datagrok-url>'
  EXTERNAL_OAUTH_RSA_PUBLIC_KEY = '<public-key-single-line-no-headers>'
  EXTERNAL_OAUTH_AUDIENCE_LIST = ('<audience-value>')
  EXTERNAL_OAUTH_TOKEN_USER_MAPPING_CLAIM = 'email'
  EXTERNAL_OAUTH_SNOWFLAKE_USER_MAPPING_ATTRIBUTE = 'LOGIN_NAME'
  EXTERNAL_OAUTH_SCOPE_MAPPING_ATTRIBUTE = 'scp'
  EXTERNAL_OAUTH_ANY_ROLE_MODE = 'ENABLE';
```

| Parameter                       | Value                                                                                          |
|---------------------------------|------------------------------------------------------------------------------------------------|
| `EXTERNAL_OAUTH_ISSUER`         | Your Datagrok instance URL, e.g. `https://datagrok.company.com`                                |
| `EXTERNAL_OAUTH_RSA_PUBLIC_KEY` | The public key from Step 1, on a single line, without `-----BEGIN/END PUBLIC KEY-----` headers |
| `EXTERNAL_OAUTH_AUDIENCE_LIST`  | Must match the "JWT Audience" field in the Datagrok connection                                 |

### Step 3: Create Snowflake Users

Each Datagrok user who will access Snowflake must have a corresponding Snowflake user with `LOGIN_NAME` matching their SSO email:

```sql
CREATE USER john
  LOGIN_NAME = 'john@company.com'
  EMAIL = 'john@company.com'
  DEFAULT_ROLE = 'ANALYST';

GRANT ROLE ANALYST TO USER john;
```

### Step 4: Configure Datagrok Connection

1. In Datagrok, create or edit a Snowflake connection
2. Fill in **Account Locator**, **Database**, **Warehouse**, **Role** as usual
3. Under authentication, select **OAuth (JWT)**
4. Upload the private key file (`datagrok_snowflake.pem`) from Step 1
5. Set **JWT Audience** to the same value used in `EXTERNAL_OAUTH_AUDIENCE_LIST`
6. Save the connection

### Step 5: Verify

Run the following query through the Datagrok connection:

```sql
SELECT CURRENT_USER(), CURRENT_ROLE()
```

It should return the Snowflake user corresponding to the logged-in Datagrok user's email.

## How Row-Level Security Works

Once OAuth (JWT) is configured, Snowflake sees each Datagrok user as a distinct identity. Row access policies that reference `CURRENT_USER()` or `CURRENT_ROLE()` apply automatically. For example:

```sql
-- Create a row access policy
CREATE ROW ACCESS POLICY user_data_policy AS (email VARCHAR) RETURNS BOOLEAN ->
  email = CURRENT_USER();

-- Apply it to a table
ALTER TABLE sensitive_data ADD ROW ACCESS POLICY user_data_policy ON (user_email);
```

Now when John queries `sensitive_data` through Datagrok, he only sees rows where `user_email` matches his identity.

## Troubleshooting

| Error                                                    | Cause                                                | Fix                                                                                       |
|----------------------------------------------------------|------------------------------------------------------|-------------------------------------------------------------------------------------------|
| `Invalid OAuth access token`                             | JWT signature validation failed                      | Verify the public key in Snowflake matches the private key uploaded to Datagrok           |
| `Incorrect username or password`                         | No Snowflake user matching the email                 | Create a Snowflake user with `LOGIN_NAME` = the user's SSO email                          |
| `The role requested is not listed in the Access Token`   | Role/scope issue                                     | Verify `EXTERNAL_OAUTH_ANY_ROLE_MODE = 'ENABLE'` and that the role is granted to the user |
| `Private key is required for OAuth (JWT) authentication` | No key uploaded                                      | Upload the private key PEM in the Datagrok connection credentials                         |
| JWT issuer mismatch                                      | `EXTERNAL_OAUTH_ISSUER` doesn't match Datagrok's URL | Set `EXTERNAL_OAUTH_ISSUER` to your exact Datagrok URL (including protocol and port)      |

## Security Notes

- The private key should be treated as a secret — anyone with it can mint tokens for any user
- JWTs are short-lived (10 minutes) and not reused
- The private key is stored encrypted in Datagrok's credentials storage
