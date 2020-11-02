import boto3
import csv
import datetime
import io
import json
import requests


# NLP package settings
region = 'us-east-1'
user = 'nlp-user'
bucket = 'nlp-bucket'
policy = 'nlp-user-policy'

# Creation of a low-level client requires an administrator account configuration.
# Run `aws configure` beforehand if you have aws-cli, or pass the additional
# parameters `aws_access_key_id` and `aws_secret_access_key` in the next line.
iam = boto3.client('iam')

# Create an S3 bucket in one of the available regions
# https://docs.aws.amazon.com/general/latest/gr/rande.html#s3_region
s3 = boto3.client('s3', region_name=region)
location = {'LocationConstraint': region}
bucket_response = s3.create_bucket(Bucket=bucket, CreateBucketConfiguration=location)

# Create an IAM user and attach a custom policy to it
user_response = iam.create_user(UserName=user)
policy_json = open('nlp-user-policy.json')
policy_document = json.load(policy_json)
policy_json.close()
policy_response = iam.create_policy(PolicyName=policy,
                                    PolicyDocument=json.dumps(policy_document),
                                    Description='Grants full access to AWS Translate and \
                                    Comprehend Medical services, and one S3 output bucket')
policy_arn = policy_response['Policy']['Arn']
iam.attach_user_policy(UserName=user, PolicyArn=policy_arn)

# Set a new access key
key_response = iam.create_access_key(UserName=user)
credentials = {
    'accessKeyId': key_response['AccessKey']['AccessKeyId'],
    'secretAccessKey': key_response['AccessKey']['SecretAccessKey']
}

# Push the active key to Datagrok's credentials storage
# https://datagrok.ai/help/govern/security#credentials-storage
grok_host = 'https://dev.datagrok.ai'
package_name = 'Nlp'
url = f'{grok_host}/api/credentials/for/{package_name}'
# Get an API Token for a special service user via `Manage | Users | Add Service User`
api_key = ''
headers = {'Authorization': api_key, 'Content-Type': 'application/json'}
cred_response = requests.post(url, json=credentials, headers=headers)
print('The credentials were added successfully' if cred_response else cred_response.status_code)

# Regularly update the access key:
# Find out when the access key was rotated last time
report_response = iam.generate_credential_report()
report = iam.get_credential_report()
table = csv.DictReader(io.StringIO(report['Content'].decode()))
for row in table:
    if row['user'] == user:
        last_rotated_iso_date = row['access_key_1_last_rotated']
        break
# Set the key activity limit
n_days = 14
last_rotated_date = datetime.datetime.fromisoformat(last_rotated_iso_date)
current_date = datetime.datetime.now(datetime.timezone.utc)
# Change the access key
if n_days < (current_date - last_rotated_date).days:
    iam.update_access_key(UserName=user,
                          AccessKeyId=credentials['accessKeyId'],
                          Status='Inactive')
    key_response = iam.create_access_key(UserName=user)
    credentials = {
        'accessKeyId': key_response['AccessKey']['AccessKeyId'],
        'secretAccessKey': key_response['AccessKey']['SecretAccessKey']
    }

# Delete the user and all related items
# https://docs.aws.amazon.com/IAM/latest/UserGuide/id_users_manage.html#id_users_deleting_cli
iam.detach_user_policy(UserName=user, PolicyArn=policy_arn)
iam.delete_policy(PolicyArn=policy_arn)
user_keys = iam.list_access_keys(UserName=user)['AccessKeyMetadata']
for key in user_keys:
    iam.delete_access_key(UserName=user, AccessKeyId=key['AccessKeyId'])
iam.delete_user(UserName=user)
