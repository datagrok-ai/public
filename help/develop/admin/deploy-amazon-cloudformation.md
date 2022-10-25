<!-- TITLE: Deployment on AWS ECS using CloudFormation -->
<!-- SUBTITLE: -->

# Deployment on AWS ECS using CloudFormation

Datagrok consist of Docker containers, [database](infrastructure.md#database)
and [persistent file storage](infrastructure.md#storage).

This document contains instructions to deploy Datagrok using [CloudFormation](https://aws.amazon.com/cloudformation/)
on [AWS ECS cluster](https://aws.amazon.com/ecs/) with [AWS RDS](https://aws.amazon.com/rds/)
and [AWS S3](https://aws.amazon.com/s3/).

We considered a lot of typical security nuances during the CloudFormation template development. As a result, you will
create a Datagrok infrastructure in AWS which applies to all standard security policies.

More information about Datagrok design and components:

* [Architecture](architecture.md)
* [Infrastructure](infrastructure.md)

## Prerequisites

1. Check that you
   have [required permissions](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/iam.list)
   on AWS account to perform CloudFormation deployment to ECS.

2. Create secret in [AWS Secret Manager](https://aws.amazon.com/secrets-manager/) with Docker Hub credentials
    1. [Create access token in Docker Hub](https://docs.docker.com/docker-hub/access-tokens/)
    2. [Create secret in Secret Manager](https://docs.aws.amazon.com/secretsmanager/latest/userguide/manage_create-basic-secret.html)
       with Docker Hub username as 'Username' and Access Token as 'Password'
    3. Copy AWS ARN for created secret. It should look like
       this: `arn:aws:secretsmanager:<region>:<account-id>:secret:<secret_name>-<random_id>`.
        * To get ARN from command line
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
           which will be valid for both endpoint: `DATAGROK_DNS`, `CVM_DNS`.
        2. Copy AWS ARN for created certificate. It should look like
           this: `arn:aws:acm:<region>:<account_id>:certificate/<certificate_id>`.
    * If you do not use AWS ACM service for SSL certificates, you can create certificate(s) for `DATAGROK_DNS`
      , `CVM_DNS` endpoints any way you are already using.
        1. [Upload certificate(s) to AWS ACM](https://docs.aws.amazon.com/acm/latest/userguide/import-certificate-api-cli.html)
        2. Copy AWS ARN for created certificate(s). It should look like
           this: `arn:aws:acm:<region>:<account_id>:certificate/<certificate_id>`.

## Deploy Datagrok components

1. Download CloudFormation
   Template: [link](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/cloudformation/cloudformation.json)
   .

2. Create CloudFormation stack
   using [AWS Console](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-console-create-stack.html)
   or [AWS CLI](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/using-cfn-cli-creating-stack.html)
    1. Use CloudFormation Template downloaded on the first step
       as [stack template](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-using-console-create-stack-template.html)

    2. [Specify stack name](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-using-console-create-stack-parameters.html)
       . It must be shorter than ten symbols to meet AWS naming requirements

    3. [Specify parameters](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-using-console-create-stack-parameters.html)
       for the stack:
        1. `ArnCvmCertificate`: Specify AWS ACM ARN for `CVM_DNS` from the 4th prerequisites step. It can be the same as
           for `ArnDatagrokCertificate`.

        2. `ArnDatagrokCertificate`: Specify AWS ACM ARN for `DATAGROK_DNS` from the 4th prerequisites step. It can be
           the same as for `ArnCvmCertificate`.

        3. `ArnDockerHubCredential`: Specify AWS Secret Manager ARN for Docker Hub authorization from 3rd prerequisites
           step.

        4. `CreateDemoData`: Datagrok provides demo databases with demo data for the full experience. Choose `true` to
           create demo databases near Datagrok.

        5. All other parameters are for Datagrok Docker images tags. The default value is `latest`.
            1. [DatagrokVersion](https://hub.docker.com/r/datagrok/datagrok)
            2. [GrokComputeVersion](https://hub.docker.com/r/datagrok/grok_connect)
            3. [H2oVersion](https://hub.docker.com/r/datagrok/h2o)
            4. [JupyterKernelGatewayVersion](https://hub.docker.com/r/datagrok/jupyter_kernel_gateway)
            5. [JupyterNotebookVersion](https://hub.docker.com/r/datagrok/jupyter_notebook)
    4. You can skip stack options; the default values suit the needs.

3. CloudFormation Stack creation takes around 10 minutes. It will create RDS, S3, ECS Cluster, and other required
   resources.

4. After the Datagrok container starts, the Datagrok server will deploy the database. You can check the status by
   checking running task log in [CloudWatch](https://aws.amazon.com/cloudwatch/)

5. Create `DATAGROK_DNS` and `CVM_DNS` DNS records which will route to the newly created Internet-facing Application
   Load Balancers.

## Configure Datagrok settings

1. Go in the web browser to `DATAGROK_DNS`, login to Datagrok using username `admin` and password `admin`.
2. Edit settings in the Datagrok (Tools | Settings...). Do not forget to click Apply to save new settings.

* Scripting:
  * CVM URL Client: `https://<CVM_DNS>`
* Dev:
  * CVM Url: `https://<CVM_DNS>`
