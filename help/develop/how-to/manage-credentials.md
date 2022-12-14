<!-- TITLE: Manage credentials -->

# Working with credentials

[Security credentials](../../govern/security.md#credentials) are used to gain access to external resources. For example,
database connections typically require a pair of login and password, whereas web services normally expect a token or an
access key. You can associate such data with two types of [entities](../../datagrok/objects.md) within the
platform: [packages](../develop.md#packages)
and [connections](../../access/data-connection.md). This article describes how to transfer this information to the
platform and common practices to keep it secure.

## Package credentials

You can store text information bound to a package in a secure manner. The easiest way to do that comes with the user
interface. Right-click on a package in Datagrok's [Packages](https://public.datagrok.ai/packages) view and select
the `Credentials...`
command from the package context menu. This will open a dialog where you can add credentials as key-value pairs. Pay
attention to the `Credentials owner` field: it may include a user or a user group, such as the current user or all users
respectively. Once added, these key-value pairs can only be read by members of the owner group.

To make a resource available to all users (the default group) running your platform extension, you can programmatically
set the relevant credentials. This is done by sending a POST request
to `${GROK_HOST}/api/credentials/for/${PACKAGE_NAME}`. The headers should contain an API key for authorization (
available on your profile page, e.g., [https://public.datagrok.ai/u](https://public.datagrok.ai/u)), and the type of
data in the request body should be `json`. For example, the following code pushes a login/password pair to
Datagrok's [credentials storage](https://datagrok.ai/help/govern/security#credentials-storage):

```javascript
let credentials = { login: 'login', password: 'password' };
let apiKey = '';  // shouldn't be empty

let grokHost = 'https://public.datagrok.ai';
let packageName = '';  // shouldn't be empty
let url = `${grokHost}/api/credentials/for/${packageName}`;

let headers = { 'Authorization': apiKey, 'Content-Type': 'application/json' };

fetch(url, { method: 'POST', headers: headers, body: JSON.stringify(credentials) })
  .then(response => console.log(response.ok));
```

Sometimes it might be more convenient to have a special service user and provide its API key in further requests. To
create one, open `Manage | Users | Actions | Add Service User`. There you can specify a login for the service user and
generate a new key.

Reading credentials is just as important. If you need them, there is a way to obtain a credentials object in your
package code:

```javascript
let _package = new DG.Package();

async function getCredentials() {
    let credentialsResponse = await _package.getCredentials();
    if (credentialsResponse === null) {
        grok.shell.info("Credentials are not set.");
        return {};
    }
    return credentialsResponse.parameters;
}
```

Additionally, you can retrieve the value of a particular parameter like this:

```javascript
_package.getCredentials().then(c => grok.shell.info(c.parameters['test']));
```

Check out this example in our [API samples](https://public.datagrok.ai/js/samples/misc/package-credentials). And to see
the full cycle of adding and reading credentials, have a look at the
public [NLP](https://github.com/datagrok-ai/public/tree/master/packages/NLP) package. It illustrates how to set new
credentials, e.g. access keys, in
a [Python script](https://github.com/datagrok-ai/public/blob/master/packages/NLP/aws/nlp-user.py)
and reach them later from the [main file](https://github.com/datagrok-ai/public/blob/master/packages/NLP/src/package.ts)
of the package.

## Database connection credentials

There are many ways to specify credentials for a database connection, namely:

* Users can go to `Data | Databases` and right-click on the data source of interest to add a connection of the specific
  provider. More generally, a data connection can be created from `Actions | Add New Connection`. Credentials can be
  passed along with other parameters in the UI.
* An equivalent way to establish a connection is to configure parameters in a `json` file. It contains a special field:

  ```json
  "credentials": {
    "parameters": {
      "login": "login",
      "password": "password"
    }
  }
  ```

  What is less noticeable from the interface, and becomes apparent here, is that credentials as parameters form a
  distinct group. You can store them in this file, the only shortcoming appears when you want to share it with others,
  e.g., push it to the repository: you have to make sure first that only your team can see these credentials there and
  not someone else who is not supposed to.
* A connection string, which lists all the parameters required for the connection and can be provided both from the UI
  and in `json`, might include a login and password, since its content is not limited in any way. However, this is not
  recommended.
* Lastly, there is an endpoint for connections that are part of a
  package: `${GROK_HOST}/api/credentials/for/${PACKAGE_NAME}.${CONNECTION_NAME}`. The rules are the same as for package
  credentials, see the example request above. This is one of the most reliable and safest options.

See also:

* [Security](../../govern/security.md#credentials)
* [How to access data](access-data.md)
* [Packages](../develop.md)
* [Data connection](../../access/data-connection.md)
