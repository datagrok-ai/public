#!/bin/bash

set -e
#set -x

function input() {
  # Check users input
  # '<secret>' '<description>' '<question>' '<options>' '<extra_values>'
  local secret=$1
  local description=$2
  local question=$3
  local options=$4
  local extra_values=$5

  [[ $extra_values == '' ]] && extra_values='not specified'

  local result=''
  [[ $description == '' ]] || echo -e "$description" >&2

  while [[ ", ${options,,}, " != *", ${result,,}, "* ]] && [[ ", ${extra_values,,}, " != *", ${result,,}, "* ]] || [[ $result == '' ]]; do
    [[ $result == '' ]] || echo "Entered value '$result' is not supported" >&2
    [[ $options == '' ]] || echo "Insert one of the following options: $options" >&2
    if [[ $secret == 'secret' ]]; then
      read -s -p "$question: " -r result
    else
      read -p "$question: " -r result
    fi
    result=$(tr '[:upper:]' '[:lower:]' <<<"$result")
    [[ $options == '' ]] && options=$result
  done
  echo "$result"
}

function check_aws_cli {
  until aws --version >/dev/null 2>&1; do
    echo 'AWS CLI is required to be installed'
    echo 'Install AWS CLI following the official documentation'
    echo 'https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html'
    echo 'Script will wait for 1 minute and check the installation again...'
    sleep 60
  done
}

echo 'Welcome to the Datagrok Deploy script!'
echo 'It will guide you through the process of deploy and create a configuration based on your needs.'
echo 'Answer the script questions and it will do the rest.'

echo -e "\n" >&2
location=$(input '' "Where would you like to setup the Datagrok application?" \
  'Datagrok location' \
  'Virtual Machine, Kubernetes, ECS' \
  'vm, k8s')

case $location in
"kubernetes" | "k8s")
  k8s_enabled='true'
  ;;
"virtual machine" | "vm")
  vm_enabled='true'
  ;;
"ecs")
  ecs_enabled='true'
  ;;
*)
  echo "Something went wrong. Contact Datagrok support."
  echo "Entered location '$location' is not supported"
  exit 1
  ;;
esac

if [[ $k8s_enabled == 'true' ]]; then
  echo 'Kubernetes is not supported by the deploy script yet.'
  echo 'Use the manual installation or reach out Datagrok support for more information.'
  exit 0
fi

if [[ $ecs_enabled == 'true' ]]; then

  echo -e "\n" >&2
  #  r53=$(input '' "Do you use Route53 for DNS management?" \
  #    'Route53 usage' \
  #    'Yes, No')
  database='rds'
  storage='s3'
  iam='role'

  if ! docker context list | grep -q datagrokAWS; then
    docker context create ecs --from-env datagrokAWS
  fi

  curl -s https://raw.githubusercontent.com/datagrok-ai/public/master/docker/ecs.datagrok.docker-compose.yaml \
    -o "datagrok.${location}.docker-compose.yaml"
  curl -s https://raw.githubusercontent.com/datagrok-ai/public/master/docker/ecs.cvm.docker-compose.yaml \
    -o "datagrok.cvm.${location}.docker-compose.yaml"
fi

if [[ $vm_enabled == 'true' ]]; then
  echo -e "\n" >&2
  setup=$(input '' "Would you like to setup Datagrok and Compute components on the same machine?\nWe recommend separate machines for Datagrok and Compute components." \
    'Datagrok setup' \
    'same, separate')

  echo -e "\n" >&2
  datagrok_host=$(input '' "" "Enter Datagrok host IP or hostname" "")

  echo -e "\n" >&2
  confirm=$(input '' "To deploy application it is required:\n1) To have access to the Datagrok host ($datagrok_host) through SSH using SSH keys\n2) Your user should be in 'docker' group on this host" \
    "Do you meet the requirements? Enter 'Confirm'" \
    'confirm')
  confirm=''

  if ! docker context list | grep -q datagrok; then
    docker context create --docker "host=ssh://${datagrok_host}:22" datagrok
  fi

  if [[ $setup == 'separate' ]]; then
    echo -e "\n" >&2
    read -p "Enter Compute host IP or hostname: " -r cvm_host

    echo -e "\n" >&2
    confirm=$(input '' "To deploy application it is required:\n1) To have access to the Datagrok host ($cvm_host) through SSH using SSH keys\n2) Your user should be in 'docker' group on this host" \
      "Do you meet the requirements? Enter 'Confirm'" \
      'confirm')
    confirm=''

    if ! docker context list | grep -q cvm; then
      docker context create --docker "host=ssh://${cvm_host}:22" cvm
    fi
  fi

  curl -s https://raw.githubusercontent.com/datagrok-ai/public/master/docker/localhost.docker-compose.yaml \
    -o "datagrok.${location}.docker-compose.yaml"
fi

echo -e "\n" >&2
adminPassword=$(input 'secret' "" "Enter Datagrok admin user password to access platform. It can not be empty" "")

db='datagrok'

if [[ $database == '' ]]; then
  echo -e "\n" >&2
  database=$(input '' 'What would you like to use as PostgreSQL database?\nDatagrok supports any PostgreSQL cluster, including cloud solutions, for example AWS RDS' \
    'PostgreSQL database setup' \
    'RDS, other')
fi

if [[ $database == 'rds' ]]; then
  echo -e "\n" >&2
  rds_create=$(input '' 'Would you like to create new RDS for Datagrok?' \
    'Create new RDS' \
    'Yes, No')
fi

if [[ $rds_create == 'yes' ]]; then
  echo -e "\n" >&2
  echo -e "Preparing RDS for Datagrok..." >&2

  check_aws_cli

  if [[ $AWS_ACCESS_KEY_ID == '' ]] || [[ $AWS_SECRET_ACCESS_KEY == '' ]]; then
    read -p "Enter AWS access key associated with an IAM user to create S3 bucket: " -r AWS_ACCESS_KEY_ID
    read -p "Enter AWS secret key associated with the access key to create S3 bucket: " -r AWS_SECRET_ACCESS_KEY
    export AWS_ACCESS_KEY_ID
    export AWS_SECRET_ACCESS_KEY
  fi

  echo -e "\n" >&2
  if [[ $AWS_REGION == '' ]]; then
    amazonRDSRegion=$(input '' "" "Enter RDS region" "")
  else
    amazonRDSRegion=$AWS_REGION
  fi
  amazonRDSName=$(input '' "" "Enter RDS name" "")

  dbAdminLogin=$(input '' "" "Enter RDS admin username" "")
  dbAdminPassword=$(input 'secret' "" "Enter RDS admin password" "")

  dbPort=5432

  echo -e "\n" >&2
  echo -e "Creating RDS for Datagrok '$amazonRDSName'..." >&2
  aws rds create-db-instance \
    --db-instance-identifier "$amazonRDSName" \
    --db-name "$db" \
    --engine 'postgres' \
    --engine-version '12.8' \
    --auto-minor-version-upgrade \
    --allocated-storage 50 \
    --max-allocated-storage 100 \
    --db-instance-class 'db.t3.large' \
    --master-username "$dbAdminLogin" \
    --master-user-password "$dbAdminPassword" \
    --port "$dbPort" \
    --multi-az \
    --no-publicly-accessible \
    --storage-encrypted \
    --deletion-protection \
    --backup-retention-period 3 \
    --tags Key=project,Value=datagrok \
    --region "${amazonRDSRegion}" \
    --output text --query 'DBInstance.[DBInstanceIdentifier, DBInstanceStatus]' || {
    rds=$(aws rds describe-db-instances --db-instance-identifier "$amazonRDSName" \
      --output text --query 'DBInstances[].[DBInstanceStatus, Endpoint.Address]')
    #    if [[ $rds == '' ]]; then
    #      echo 'Something went wrong during RDS creation. Check the log above for an error.' >&2
    #      exit 1
    #    fi
    rds_status=$(awk '{print $1}' <<<"$rds")
    echo -e "\n" >&2
    echo "RDS '$amazonRDSName' already exists" >&2
  }
  until [[ $rds_status == 'available' ]]; do
    echo "Waiting 1 minute for RDS for Datagrok '$amazonRDSName' to become available..." >&2
    rds=$(aws rds describe-db-instances --db-instance-identifier "$amazonRDSName" \
      --output text --query 'DBInstances[].[DBInstanceStatus, Endpoint.Address]')
    rds_status=$(awk '{print $1}' <<<"$rds")
    sleep 60
  done

  dbServer=$(awk '{print $2}' <<<"$rds")
else
  echo -e "\n" >&2
  echo 'Provide Datagrok internal PostgreSQL credentials' >&2
  dbServer=$(input '' "" "Enter Datagrok internal database server address" "")
  dbPort=$(input '' "" "Enter Datagrok internal database server port" "")

  echo -e "\n" >&2
  echo 'Provide Datagrok admin credentials to connect to internal database' >&2
  echo 'It will be used once to create required schema for Datagrok application' >&2
  dbAdminLogin=$(input '' "" "Enter Postgres admin username" "")
  dbAdminPassword=$(input 'secret' "" "Enter Postgres admin password" "")
fi

#echo -e "\n" >&2
#confirm=$(input "Is the connection to the internal Datagrok database TLS encrypted?" \
#  "TLS encryption" \
#  'yes, no')
#dbSsl=$confirm
#confirm=''

echo -e "\n" >&2
echo 'Choose Datagrok credentials to connect to internal database' >&2
dbLogin=$(input '' "" "Enter Datagrok username to connect to internal database" "")
dbPassword=$(input 'secret' "" "Enter Datagrok password to connect to internal database" "")

GROK_PARAMETERS="{\\\\\"dbServer\\\\\": \\\\\"${dbServer}\\\\\",\\\\\"db\\\\\": \\\\\"${db}\\\\\",\\\\\"dbAdminLogin\\\\\":\\\\\"${dbAdminLogin}\\\\\",\\\\\"dbAdminPassword\\\\\": \\\\\"${dbAdminPassword}\\\\\",\\\\\"dbLogin\\\\\": \\\\\"${dbLogin}\\\\\",\\\\\"dbPassword\\\\\": \\\\\"${dbPassword}\\\\\",\\\\\"adminPassword\\\\\": \\\\\"${adminPassword}\\\\\""

if [[ $storage == '' ]]; then
  echo -e "\n" >&2
  storage=$(input '' 'What would you like to use as persistent storage?' \
    'Persistent storage location' \
    'Local File System, S3, GCS' \
    'local, fs, local fs, system, file system')
fi

if [[ $storage == 's3' ]]; then
  echo -e "\n" >&2
  s3_create=$(input '' 'Would you like to create new S3 bucket for Datagrok?' \
    'Create new S3 bucket' \
    'Yes, No')

  echo -e "\n" >&2
  if [[ $s3_create == 'no' ]]; then
    amazonStorageRegion=$(input '' "" "Enter S3 region" "")
  elif [[ $AWS_REGION == '' ]]; then
    amazonStorageRegion=$(input '' "" "Enter S3 region" "")
  else
    amazonStorageRegion=$AWS_REGION
  fi

  amazonStorageBucket=$(input '' "" "Enter S3 bucket name" "")
  GROK_PARAMETERS+=",\\\\\"amazonStorageRegion\\\\\": \\\\\"${amazonStorageRegion}\\\\\",\\\\\"amazonStorageBucket\\\\\": \\\\\"${amazonStorageBucket}\\\\\""

  if [[ $s3_create == 'yes' ]]; then
    echo -e "\n" >&2
    echo -e "Creating S3 bucket for Datagrok..." >&2

    check_aws_cli

    if [[ $AWS_ACCESS_KEY_ID == '' ]] || [[ $AWS_SECRET_ACCESS_KEY == '' ]]; then
      read -p "Enter AWS access key associated with an IAM user to create S3 bucket: " -r AWS_ACCESS_KEY_ID
      read -p "Enter AWS secret key associated with the access key to create S3 bucket: " -r AWS_SECRET_ACCESS_KEY
      export AWS_ACCESS_KEY_ID
      export AWS_SECRET_ACCESS_KEY
    fi

    if aws s3api head-bucket --bucket "$amazonStorageBucket" 2>/dev/null; then
      echo -e "\n" >&2
      echo "'$amazonStorageBucket' S3 bucket already exists" >&2
    else
      aws s3api create-bucket \
        --bucket "$amazonStorageBucket" \
        --create-bucket-configuration LocationConstraint="$amazonStorageRegion" \
        --region "$amazonStorageRegion"
    fi
    aws s3api put-bucket-encryption \
      --bucket "$amazonStorageBucket" \
      --server-side-encryption-configuration '{"Rules": [{"ApplyServerSideEncryptionByDefault": {"SSEAlgorithm": "AES256"}}]}'
    aws s3api put-public-access-block \
      --bucket "$amazonStorageBucket" \
      --public-access-block-configuration "BlockPublicAcls=true,IgnorePublicAcls=true,BlockPublicPolicy=true,RestrictPublicBuckets=true"
  fi

  if [[ $iam == '' ]]; then
    echo -e "\n" >&2
    iam=$(input '' 'Would you like to use IAM credentials or IAM role for Datagrok to access S3 bucket?\nWe recommend to use IAM roles on AWS resources if possible.' \
      'IAM role or IAM credentials' \
      'IAM Role, IAM Credentials' \
      'role, credentials')
  fi

  if [[ $iam == 'credentials' ]]; then
    read -p "Enter S3 credential ID: " -r amazonStorageId
    read -p "Enter S3 credential secret key: " -r amazonStorageKey

    GROK_PARAMETERS+=",\\\\\"amazonStorageId\\\\\": \\\\\"${amazonStorageId}\\\\\",\\\\\"amazonStorageKey\\\\\": \\\\\"${amazonStorageKey}\\\\\""
  fi
fi

if [[ $storage == 'gcs' ]]; then
  echo -e "\n" >&2
  googleStorageCert=$(input '' "" "Enter Access certificate to Google Cloud Storage" "")

  GROK_PARAMETERS+=",\\\\\"googleStorageCert\\\\\": \\\\\"${googleStorageCert}\\\\\""
fi

echo
echo 'Creating Datagrok configuration from input data...'
echo
GROK_PARAMETERS+="}"

perl -i -p0e "s%GROK_PARAMETERS: \"{.*}\"%GROK_PARAMETERS: \"${GROK_PARAMETERS}\"%s" "datagrok.${location}.docker-compose.yaml"

echo
echo "Deploying Datagrok to ${location}..."
echo

case $location in
"kubernetes" | "k8s")
  echo 'K8S not supported yet'
  ;;

"virtual machine" | "vm")
  if [[ $setup == 'separate' ]]; then
    docker context use datagrok >/dev/null 2>&1
    docker-compose --project-name datagrok \
      --profile datagrok \
      -f "datagrok.${location}.docker-compose.yaml" up -d
    docker context use cvm >/dev/null 2>&1
    docker-compose --project-name cvm \
      --profile cvm \
      -f "datagrok.${location}.docker-compose.yaml" up -d
  else
    docker context use datagrok >/dev/null 2>&1
    docker-compose \
      --project-name datagrok \
      --profile datagrok --profile cvm \
      -f "datagrok.${location}.docker-compose.yaml" up -d
  fi
  docker context use default >/dev/null 2>&1
  rm "datagrok.${location}.docker-compose.yaml"
  ;;

"ecs")
  docker context use datagrokAWS >/dev/null 2>&1
  docker compose --project-name datagrok -f "datagrok.${location}.docker-compose.yaml" up || true
  if [[ $rds_create == 'yes' ]]; then
    ecs_sg=$(aws ecs describe-services \
      --cluster datagrok \
      --services "$(aws ecs list-services --cluster datagrok --output text --query 'serviceArns[0]')" \
      --output text --query 'services[].networkConfiguration.awsvpcConfiguration.securityGroups[]')
    aws rds modify-db-instance \
      --db-instance-identifier "${amazonRDSName}" \
      --region "${amazonRDSRegion}" \
      --vpc-security-group-ids "${ecs_sg}" \
      --apply-immediately --output text --query '[DBInstance.DBInstanceIdentifier, VpcSecurityGroups[]]'
  fi
  docker compose --project-name datagrok-cvm -f "datagrok.cvm.${location}.docker-compose.yaml" up || true
  docker context use default >/dev/null 2>&1
  rm "datagrok.${location}.docker-compose.yaml" "datagrok.cvm.${location}.docker-compose.yaml"
  ;;

*)
  echo "Something went wrong. Contact Datagrok support."
  echo "Entered location '$location' is not supported"
  exit 1
  ;;
esac
