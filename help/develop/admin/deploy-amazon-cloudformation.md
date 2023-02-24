---
title: "Deployment on AWS ECS using CloudFormation"
---

Datagrok is based on Docker containers, [database](infrastructure.md#database)
and [persistent file storage](infrastructure.md#storage).

This document contains instructions to deploy Datagrok using [CloudFormation](https://aws.amazon.com/cloudformation/)
on [AWS ECS cluster](https://aws.amazon.com/ecs/) with [AWS RDS](https://aws.amazon.com/rds/)
and [AWS S3](https://aws.amazon.com/s3/).

We considered a lot of typical security nuances during the CloudFormation template development. As a result, you will
create a Datagrok infrastructure in AWS that applies to all standard security policies.

More information about Datagrok design and components:

* [Architecture](architecture.md)
* [Infrastructure](infrastructure.md)

## Prerequisites

1. Check that you
   have [required permissions](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/iam.list)
   on AWS account to perform CloudFormation deployment to ECS.

2. Create a secret in [AWS Secret Manager](https://aws.amazon.com/secrets-manager/) with Docker Hub credentials
   1. [Create access token in Docker Hub](https://docs.docker.com/docker-hub/access-tokens/)
   2. [Create a secret in Secret Manager](https://docs.aws.amazon.com/secretsmanager/latest/userguide/manage_create-basic-secret.html)
      with Docker Hub username as 'Username' and Access Token as 'Password'
   3. Copy AWS ARN for the created secret. It should look like this:
      `arn:aws:secretsmanager:<region>:<account-id>:secret:<secret_name>-<random_id>`.

        * To get ARN from the command line
          using [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) run:

          ```shell
          aws secretsmanager list-secrets --filters Key=name,Values=<secret_name> --output text --query 'SecretList[].ARN'
          ```

3. Come up with two endpoints: `DATAGROK_DNS`, `CVM_DNS`. Datagrok requires two endpoints: `DATAGROK_DNS` and `CVM_DNS`.
   Users will use `DATAGROK_DNS` to access Datagrok Web UI, and requests `CVM_DNS` will be sent automatically by
   Datagrok Client.

4. Create RSA SSL certificate(s) for `DATAGROK_DNS`, `CVM_DNS`.

   * If you use AWS ACM service for SSL certificates
       1. [Generate ACM certificate in AWS](https://docs.aws.amazon.com/acm/latest/userguide/gs-acm-request-public.html)
          which will be valid for both endpoints: `DATAGROK_DNS`, `CVM_DNS`.
       2. Copy AWS ARN for the created certificate. It should look like
          this: `arn:aws:acm:<region>:<account_id>:certificate/<certificate_id>`.
   * If you do not use AWS ACM service for SSL certificates, you can create a certificate(s) for `DATAGROK_DNS`
     , `CVM_DNS` endpoints any way you are already using.
       1. [Upload certificate(s) to AWS ACM](https://docs.aws.amazon.com/acm/latest/userguide/import-certificate-api-cli.html)
       2. Copy AWS ARN for the created certificate(s). It should look like
          this: `arn:aws:acm:<region>:<account_id>:certificate/<certificate_id>`.

## Deploy Datagrok components

1. Download CloudFormation
   Template in
   [YAML](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/cloudformation/cloudformation.yml)
   or
   [JSON](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/cloudformation/cloudformation.json)
   format as you prefer.

2. Create a CloudFormation stack
   using [AWS Console](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-console-create-stack.html)
   or [AWS CLI](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/using-cfn-cli-creating-stack.html)

   1. Use CloudFormation Template downloaded on the first step
      as [stack template](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-using-console-create-stack-template.html)
   2. [Specify stack name](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-using-console-create-stack-parameters.html)
      . It must be shorter than ten symbols to meet AWS naming requirements
   3. [Specify parameters](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-using-console-create-stack-parameters.html)
      for the stack:
      1. `ArnCvmCertificate`: Specify AWS ACM ARN for `CVM_DNS` from the 4th prerequisites step. It can be the same as
         for the `ArnDatagrokCertificate`.
      2. `ArnDatagrokCertificate`: Specify AWS ACM ARN for `DATAGROK_DNS` from the 4th prerequisites step. It can be
         the same as for `ArnCvmCertificate`.
      3. `ArnDockerHubCredential`: Specify AWS Secret Manager ARN for Docker Hub authorization from 3rd prerequisites
         step.
      4. `CreateDemoData`: Datagrok provides demo databases with demo data for the full experience. Choose `true` to
         create demo databases near Datagrok.
      5. `LaunchType`: It is an optional argument. The default value is `FARGATE`. It will suit best in most cases. Also,
         the template supports the `EC2` launch type, which will use EC2 instances and reduce the price of the stand.
         More information about EC2 launch type is described [below](#optional-cost-reduction-stand).
      6. `Ec2PublicKey`: It is an optional argument. It is only required for `EC2` `LaunchType`. By default, you can
         skip it. More information about EC2 launch type is described [below](#optional-cost-reduction-stand).
      7. All other parameters are for Datagrok Docker image tags. The default value is `latest`.
         1. [DatagrokVersion](https://hub.docker.com/r/datagrok/datagrok)
         2. [GrokComputeVersion](https://hub.docker.com/r/datagrok/grok_connect)
         3. [H2oVersion](https://hub.docker.com/r/datagrok/h2o)
         4. [JupyterKernelGatewayVersion](https://hub.docker.com/r/datagrok/jupyter_kernel_gateway)
         5. [JupyterNotebookVersion](https://hub.docker.com/r/datagrok/jupyter_notebook)
   4. You can skip stack options; the default values suit the needs.

3. CloudFormation Stack creation takes around 10 minutes. It will create RDS, S3, ECS Cluster, and other required
   resources.

4. After the Datagrok container starts, the Datagrok server will deploy the database. You can check the status by
   checking the running task log in [CloudWatch](https://aws.amazon.com/cloudwatch/)

5. Create `DATAGROK_DNS` and `CVM_DNS` DNS records which will route to the newly created Internet-facing Application
   Load Balancers.

## Configure Datagrok settings

1. Go to the web browser to `DATAGROK_DNS`, login to Datagrok using username `admin` and password `admin`.
2. Edit settings in the Datagrok (Tools | Settings...). Remember to click Apply to save new settings.

   * Scripting:
     * Api Url: `https://<DATAGROK_DNS>`
     * Cvm Url: `https://<CVM_DNS>`
     * H2o Url: `https://<CVM_DNS>:54321`
     * Cvm Url Client: `https://<CVM_DNS>`
   * Admin:
     * Web Root: `https://<DATAGROK_DNS>`
     * Api Root: `https://<DATAGROK_DNS>/api`

3. Reload the page and re-login. Now you are good to go.

## Optional: Cost reduction stand

AWS stack uses `FARGATE` instances for deployment by default. To reduce
infrastructure costs, you can use EC2 instances. To do so, follow the instructions above with additions below in
the [prerequisites](#ec2-prerequisites) and [parameters](#ec2-parameters).

### EC2 Prerequisites

1. Before deploying the Datagrok Stand in addition to [Prerequisites](#prerequisites), create RSA key pair. It is
   required to get access to the instances that will be created, you need to have SSH key pair: a private key and a
   public key.

     * If you already have an RSA key pair, you can use the existing one.
     * If you have a Linux-based OS or macOS, type in terminal `ssh-keygen` and hit **Enter**.
       You’ll be asked to enter a passphrase. Hit **Enter** to skip this step.
       It will create `id_rsa` and `id_rsa.pub` files in the `~/.ssh` directory.
     * If you have Windows open the Settings panel, then click Apps.
       Under the *Apps and Features* heading, click **Optional Features**.
       Scroll down the list to see if **OpenSSH Client** is listed.
       If it’s not, click the plus-sign next to **Add a feature**.
       Scroll through the list to find and select **OpenSSH Client**. Finally, click **Install**.
       Press the **Windows key**, type **cmd** under *Best Match*, and right-click **Command Prompt**.
       Click **Run as Administrator**.
       If prompted, click **Yes** in the *Do you want to allow this app to make changes to your device?* pop-up.
       In the command prompt, type `ssh-keygen` and hit **Enter**.
       You’ll be asked to enter a passphrase. Hit **Enter** to skip this step.
       By default, the system will save the keys to `C:\Users\your_username/.ssh/id_rsa`.

2. Copy the content of the public part of the key pair: `id_rsa.pub`. It will be placed in the EC2 instances
   using the `Ec2PublicKey` parameter to access machines.

### EC2 Parameters

1. Change the `LaunchType` parameter to `EC2`.

2. Paste your public key content from [EC2 prerequisites](#ec2-prerequisites) to
   the parameter `Ec2PublicKey`.
