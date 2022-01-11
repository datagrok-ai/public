#!/bin/bash

#function input() {
#  # Check users input
#  # '<description>' '<question>' variable
#  local description=$1
#  local question=$2
#
#  echo "$description"
#  while [[ $location != 'virtual machine' ]] || [[ $location != 'vm' ]] || [[ $location != 'kubernetes' ]] || [[ $location != 'k8s' ]] || [[ $location != 'ecs' ]] || [[ $location == '' ]]; do
#    [[ $location == '' ]] || echo "Entered location '$location' is not supported"
#    echo 'Insert one of the following options: Virtual Machine, Kubernetes, ECS'
#    read -p 'Enter Datagrok location: ' -r location
#    location=$(tr '[:upper:]' '[:lower:]' <<<"$location")
#  done
#}

echo 'Welcome to the Datagrok Deploy script!'
echo 'It will guide you through the process of deploy and create a configuration based on your needs.'
echo 'Answer the script questions and it will do the rest.'
echo

echo 'Where would you like to setup the Datagrok application?'
while [[ $location != 'virtual machine' ]] && [[ $location != 'vm' ]] && [[ $location != 'kubernetes' ]] && [[ $location != 'k8s' ]] && [[ $location != 'ecs' ]] || [[ $location == '' ]]; do
  [[ $location == '' ]] || echo "Entered location '$location' is not supported"
  echo 'Enter one of the following options: Virtual Machine, Kubernetes, ECS'
  read -p 'Datagrok location: ' -r location
  location=$(tr '[:upper:]' '[:lower:]' <<<"$location")
done

case $location in
"kubernetes")
  k8s_enabled='true'
  ;;
"k8s")
  k8s_enabled='true'
  ;;
"virtual machine")
  vm_enabled='true'
  ;;
"vm")
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
  echo 'ECS is not supported by the deploy script yet.'
  echo 'Use the manual installation or reach out Datagrok support for more information.'
  exit 0
fi

if [[ $vm_enabled == 'true' ]]; then
  echo 'Would you like to setup Datagrok and Compute components on the same machine?'
  echo 'We recommend separate machines for Datagrok and Compute components.'
  while [[ $setup != 'same' ]] && [[ $setup != 'separate' ]] || [[ $setup == '' ]]; do
    [[ $setup == '' ]] || echo "Entered location '$setup' is not supported"
    echo 'Enter one of the following options: same, separate'
    read -p 'Datagrok setup: ' -r setup
    setup=$(tr '[:upper:]' '[:lower:]' <<<"$setup")
  done

  read -p "Enter Datagrok host IP or hostname: " -r datagrok_host
  echo 'To deploy application it is required:'
  echo "1) To have access to the Datagrok host ($datagrok_host) through SSH using SSH keys"
  echo "2) Your user should be in 'docker' group on this host"
  while [[ $confirm != 'confirm' ]] || [[ $confirm == '' ]]; do
    [[ $confirm == '' ]] || echo "Meet the requirements and come back."
    echo 'Do you meet the requirements?'
    read -p "Enter 'Confirm': " -r confirm
    confirm=$(tr '[:upper:]' '[:lower:]' <<<"$confirm")
  done
  confirm=''
  docker context create --docker "host=ssh://${datagrok_host}:22" datagrok
  if [[ $setup == 'separate' ]]; then
    read -p "Enter Compute host IP or hostname: " -r cvm_host
    echo 'To deploy application it is required:'
    echo "1) To have access to the Compute host ($cvm_host) through SSH using SSH keys"
    echo "2) Your user should be in 'docker' group on this host"
    while [[ $confirm != 'confirm' ]] || [[ $confirm == '' ]]; do
      [[ $confirm == '' ]] || echo "Meet the requirements and come back."
      echo 'Do you meet the requirements?'
      read -p "Enter 'Confirm': " -r confirm
      confirm=$(tr '[:upper:]' '[:lower:]' <<<"$confirm")
    done
    confirm=''
    docker context create --docker "host=ssh://${cvm_host}:22" cvm
  fi

  while [[ $adminPassword == '' ]]; do
    read -p "Enter Datagrok admin user password to access platform. It can not be empty: " -r adminPassword
  done

  echo 'Provide Datagrok internal PostgreSQL credentials'
  echo 'Datagrok supports any PostgreSQL cluster, including cloud solutions, for example AWS RDS'
  read -p 'Enter Datagrok internal database server address: ' -r dbServer
  read -p 'Enter Datagrok internal database server port: ' -r dbPort
  echo 'Provide Datagrok credentials to connect to internal database'
  read -p 'Enter Datagrok username to connect to internal database: ' -r dbLogin
  read -p 'Enter Datagrok password to connect to internal database: ' -r dbPassword
  while [[ $confirm != 'yes' ]] && [[ $confirm != 'no' ]] || [[ $confirm == '' ]]; do
    [[ $confirm == '' ]] || echo "Entered answer '$confirm' is not supported"
    echo 'Is the connection to the internal Datagrok database TLS encrypted?'
    read -p "Enter 'yes' or 'no': " -r confirm
    confirm=$(tr '[:upper:]' '[:lower:]' <<<"$confirm")
  done
  dbSsl=confirm
  confirm=''
  echo 'Provide Datagrok admin credentials to connect to internal database'
  echo 'It will be used once to create required schema for Datagrok application'
  read -p 'Enter Postgres admin username: ' -r dbAdminLogin
  read -p 'Enter Postgres admin password: ' -r dbAdminPassword

  GROK_PARAMETERS="{\\\\\"dbServer\\\\\": \\\\\"${dbServer}\\\\\",\\\\\"db\\\\\": \\\\\"datagrok\\\\\",\\\\\"dbAdminLogin\\\\\":\\\\\"${dbAdminLogin}\\\\\",\\\\\"dbAdminPassword\\\\\": \\\\\"${dbAdminPassword}\\\\\",\\\\\"dbLogin\\\\\": \\\\\"${dbLogin}\\\\\",\\\\\"dbPassword\\\\\": \\\\\"${dbPassword}\\\\\",\\\\\"adminPassword\\\\\": \\\\\"${adminPassword}\\\\\""

  echo 'What would you like to use as persistent storage?'
  while [[ $storage != 'local' ]] && [[ $storage != 's3' ]] && [[ $storage != 'gcs' ]] || [[ $storage == '' ]]; do
    [[ $storage == '' ]] || echo "Entered location '$location' is not supported"
    echo 'Enter one of the following options: local, s3, gcs'
    read -p 'Persistent storage location: ' -r storage
    storage=$(tr '[:upper:]' '[:lower:]' <<<"$storage")
  done

  if [[ $storage == 's3' ]]; then
    read -p "Enter S3 region: " -r amazonStorageRegion
    read -p "Enter S3 bucket name: " -r amazonStorageBucket
    GROK_PARAMETERS+=",\\\\\"amazonStorageRegion\\\\\": \\\\\"${amazonStorageRegion}\\\\\",\\\\\"amazonStorageBucket\\\\\": \\\\\"${amazonStorageBucket}\\\\\""
    read -p "Enter S3 credential ID, Datagrok will resolve IAM role if empty: " -r amazonStorageId
    read -p "Enter S3 credential secret key, Datagrok will resolve IAM role if empty: " -r amazonStorageKey

    if [[ $amazonStorageId != '' ]] && [[ $amazonStorageKey != '' ]]; then
      GROK_PARAMETERS+=",\\\\\"amazonStorageId\\\\\": \\\\\"${amazonStorageId}\\\\\",\\\\\"amazonStorageKey\\\\\": \\\\\"${amazonStorageKey}\\\\\""
    fi
  fi

  if [[ $storage == 'gcs' ]]; then
    read -p "Enter Access certificate to Google Cloud Storage: " -r googleStorageCert

    GROK_PARAMETERS+=",\\\\\"googleStorageCert\\\\\": \\\\\"${googleStorageCert}\\\\\""
  fi

  echo
  echo 'Creating Datagrok configuration from input data...'
  echo
  GROK_PARAMETERS+="}"

  curl -O https://raw.githubusercontent.com/datagrok-ai/public/master/docker/localhost.docker-compose.yaml

  if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' -e "s%GROK_PARAMETERS: \".*\"%GROK_PARAMETERS: \"${GROK_PARAMETERS}\"%g" localhost.docker-compose.yaml
  else
    sed -i -e "s%GROK_PARAMETERS: \".*\"%GROK_PARAMETERS: \"${GROK_PARAMETERS}\"%g" localhost.docker-compose.yaml
  fi

  if [[ $setup == 'separate' ]]; then
    docker context use datagrok
    docker-compose --project-name datagrok --profile datagrok -f localhost.docker-compose.yaml up -d
    docker context use cvm
    docker-compose --project-name cvm --profile cvm up -d
  else
    docker context use datagrok
    docker-compose --project-name datagrok --profile datagrok --profile cvm -f localhost.docker-compose.yaml up -d
  fi
  docker context use default
fi
