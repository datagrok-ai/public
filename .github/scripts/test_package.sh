#!/bin/bash

set -x -o pipefail

crnt=$(pwd)

PACKAGE=$1

SKIP_CLEANUP=$2

DIR="packages/$PACKAGE"

dockers=$(find ${DIR}/dockerfiles -type f -name Dockerfile | tr '\n' ' ')

profiles='--profile datagrok --profile db'
if [ -d "${DIR}/dockerfiles" ]; then
    profiles+=' --profile grok_spawner'

    for dockerfile in $(echo -e $dockers); do
      docker_config_file="$(dirname $dockerfile)/container.json"
      cp $docker_config_file $docker_config_file.bkp
      if [ -f "${docker_config_file}" ]; then
        jq 'del(.gpu)' ${docker_config_file} > temp.json && mv temp.json ${docker_config_file}
      fi
    done
fi

if [ -d "${DIR}/connections" ] || [[ -n $(find ${DIR}/node_modules/@datagrok -type d -name connections) ]] || \
   [ -d "${DIR}/queries" ] || [[ -n $(find ${DIR}/node_modules/@datagrok -type d -name queries) ]]; then
profiles+=' --profile grok_connect'
fi

grok_deps="$(jq  -r '. | select( has("devDependencies") == true ).devDependencies | to_entries[] | .key | select(test("@datagrok/.*")?)' ${DIR}/package.json)"
if [[ "$(tr '[:upper:]' '[:lower:]' <<<${PACKAGE})" == "chem" ]] || \
   [[ "$(tr '[:upper:]' '[:lower:]' <<<${PACKAGE})" == "simpkpd" ]] || \
   [[ "$(tr '[:upper:]' '[:lower:]' <<<${PACKAGE})" == "dendrogram" ]] || \
   [[ "$(tr '[:upper:]' '[:lower:]' <<<${PACKAGE})" == "cvmtests" ]] || \
   [[ "$(tr '[:upper:]' '[:lower:]' <<<${PACKAGE})" == "samples" ]] || \
   [[ "$(tr '[:upper:]' '[:lower:]' <<<${grok_deps})" == *"chem"* ]] || \
   [[ "$(tr '[:upper:]' '[:lower:]' <<<${grok_deps})" == *"simpkpd"* ]] || \
   [[ "$(tr '[:upper:]' '[:lower:]' <<<${grok_deps})" == *"cvmtests"* ]] || \
   [[ "$(tr '[:upper:]' '[:lower:]' <<<${grok_deps})" == *"dendrogram"* ]] || \
   [[ "$(tr '[:upper:]' '[:lower:]' <<<${grok_deps})" == *"samples"* ]]; then
echo "Add Scripting as dependency for package ${PACKAGE}"
profiles+=' --profile scripting'
fi

dependencies="$(jq '(. | select( has("dependencies") == true ).dependencies) * (. | select( has("devDependencies") == true ).devDependencies)' "${DIR}/package.json")"
unpublished_deps="$(jq -r '. | to_entries | map(select(.value | match("\\.\\./.*")))[] | "\(.key)=\(.value)"' <<<$dependencies | tr '\n' ' ')"
if [[ "${unpublished_deps}" == "" ]] && [[ "$(tr '[:upper:]' '[:lower:]' <<<${PACKAGE})" != *"tests"* ]]; then
    DATAGROK_VERSION='latest'
else
    DATAGROK_VERSION='bleeding-edge'
    export GROK_SPAWNER_VERSION='bleeding-edge'
    export GROK_JUPYTER_KERNEL_GATEWAY_VERSION='bleeding-edge'
fi
export DATAGROK_VERSION

export DATAGROK_PORT=8765
export GROK_SPAWNER_PORT=8785
export GROK_CONNECT_PORT=2345
export DATAGROK_DB_PORT=15434

cp docker/localhost.docker-compose.yaml docker/github_actions.docker-compose.yaml
sed -i '/"dbServer": "database"/a\ \ \ \ \ \ \ \ \ \ "initialSetupCompleted": true,' docker/github_actions.docker-compose.yaml || exit 1

datlasUrl="http://127.0.0.1:${DATAGROK_PORT}"
apiUrl="http://127.0.0.1:${DATAGROK_PORT}/api"
alias='githubactions'

# Cleanup function (like 'finally')
cleanup() {
    if [ -z "$SKIP_CLEANUP" ]; then
        echo "Performing cleanup..."
        cd ${crnt} || exit 1

        git status
        git checkout .
        git stash pop
        if [ -d "${DIR}/dockerfiles" ]; then
            for dockerfile in $(echo -e $dockers); do
              docker_config_file="$(dirname $dockerfile)/container.json"
              mv $docker_config_file.bkp $docker_config_file
            done
        fi
        docker compose -p ${alias} -f "docker/github_actions.docker-compose.yaml" --profile all down --volumes
        docker rm -f $(docker ps --format {{.Names}} | grep ${alias}) || true
        docker compose -p ${alias} -f "docker/github_actions.docker-compose.yaml" --profile all down --volumes
        rm -rf "docker/github_actions.docker-compose.yaml"
    else
      echo "Cleanup is skipped."
    fi
}

# Trap EXIT to always run the cleanup function
trap cleanup EXIT

docker compose -p ${alias} -f "docker/github_actions.docker-compose.yaml" --profile all down --volumes || exit 1
docker rm -f $(docker ps --format {{.Names}} | grep ${alias}) || true
docker compose -p ${alias} -f "docker/github_actions.docker-compose.yaml" ${profiles} pull || exit 1
docker compose -p ${alias} -f "docker/github_actions.docker-compose.yaml" ${profiles} up -d || exit 1

cd ${DIR} || exit 1
git stash save -u "${alias}"
echo 'Removing...'
git clean -ndX
git clean -fdX
npm install
if [[ "${unpublished_deps}" != "" ]]; then
    grok link
fi
npm run build -- --mode=production || exit 1

cd ${crnt} || exit 1

until .github/scripts/check-output.sh "curl -s ${apiUrl}/info/server" '"Http Server"'
do
    sleep 1
done
token=$(curl -s -X POST ${apiUrl}/users/login/dev/admin | jq -r .token)
# Check if the token is null
if [ "$token" == "null" ] || [ -z "$token" ]; then
  echo "Error: Token is null or empty."
  exit 1
else
  echo "Token is valid: $token"
fi
long_term_token=$(curl -s -X POST "${apiUrl}/users/sessions/current/refresh" -H "Authorization: ${token}" | jq -r .token)
# Check if the token is null
if [ "$long_term_token" == "null" ] || [ -z "$long_term_token" ]; then
  echo "Error: Token is null or empty."
  exit 1
else
  echo "Token is valid: $token"
fi
curl -s -H "Authorization: $long_term_token" -H "Content-Type": "application/json" "${apiUrl}/admin/plugins/admin/settings" -X POST -d "{\"#type\":\"AdminPluginSettings\", \"agreementDate\":null, \"webRoot\": \"${datlasUrl}\", \"apiRoot\": \"${apiUrl}\"}" || exit 1

key='admin'
grok config add --default --alias ${alias} --server "${apiUrl}" --key "$key" || exit 1

cd ${DIR} || exit 1
grok_deps="$(jq  -r '. | select( has("devDependencies") == true ).devDependencies | to_entries[] | .key | select(test("@datagrok/.*")?)' package.json)"
if [ -n "$grok_deps" ]; then
for dep in $grok_deps; do
  current_dir=$(pwd)
  cd $(echo "node_modules/$dep" | tr -d '\r') || exit 1
  count=0
  retries=5
  echo "Publishing $dep to ${alias}..."
  until grok publish ${alias}; do
    exit=$?
    wait=$((2 ** count))
    count=$((count + 1))
    if [ $count -lt "$retries" ]; then
      echo "Retry $count/$retries exited $exit, retrying 'grok publish ${alias}' for $dep in $wait seconds..."
      sleep $wait
    else
      echo "Retry $count/$retries exited $exit, no more retries left for 'grok publish ${alias}' for $dep."
      exit $exit
    fi
  done
  cd $current_dir || exit 1
done
fi

count=0
retries=5
until grok publish ${alias}; do
    exit=$?
    wait=$((2 ** count))
    count=$((count + 1))
    if [ "$count" -lt "$retries" ]; then
      echo "Retry $count/$retries exited $exit, retrying 'grok publish ${alias}' in $wait seconds..."
      sleep $wait
    else
      echo "Retry $count/$retries exited $exit, no more retries left for 'grok publish ${alias}'."
      exit $exit
    fi
done

cd ${crnt} || exit 1

if [[ "${PACKAGE}" == "Chemspace" ]]; then
  curl -s -H "Authorization: $long_term_token" -H "Content-Type": "application/json" "${apiUrl}/credentials/for/Chemspace.Chemspace" -X POST -d "{\"apiKey\": \"$CHEMSPACE_APIKEY\"}"
fi

check_container=false
for dockerfile in $(echo -e ${dockers}); do
docker_config_file="$(dirname $dockerfile)/container.json"
if [ -f "${docker_config_file}" ]; then
  on_demand=$(jq '.on_demand' ${docker_config_file})
  if [ "$on_demand" == "true" ]; then
    check_container=false
  else
    check_container=true
    break
  fi
else
  check_container=true
  break
fi
done


token=$(curl -s -X POST ${apiUrl}/users/login/dev/admin | jq -r .token)
# Check if the token is null
if [ "$token" == "null" ] || [ -z "$token" ]; then
  echo "Error: Token is null or empty."
  exit 1
else
  echo "Token is valid: $token"
fi
long_term_token=$(curl -s -X POST "${apiUrl}/users/sessions/current/refresh" -H "Authorization: ${token}" | jq -r .token)
# Check if the token is null
if [ "$long_term_token" == "null" ] || [ -z "$long_term_token" ]; then
  echo "Error: Token is null or empty."
  exit 1
else
  echo "Token is valid: $token"
fi

if [[ "${check_container}" == "true" ]] ; then
    if [[ "${alias}" == "test" ]]; then
      key="$TEST_TEST_DEV_KEY"
    elif [[ "${alias}" == "dev" ]]; then
      key="$DEV_TEST_DEV_KEY"
    else
      key='admin'
    fi

    count=0
    retries=20
    until .github/scripts/check-output.sh "curl -s -H 'Authorization: $long_term_token' '${apiUrl}/docker/containers'" "${PACKAGE}"
    do
      wait=$((2 ** count))
      count=$((count + 1))
      if [ "$count" -lt "$retries" ]; then
          docker compose -p ${alias} -f "docker/github_actions.docker-compose.yaml" --profile all exec db sh -c "PAGER=cat psql -U postgres -d datagrok -c 'SELECT * FROM docker_containers;'"
          sleep $wait
          echo -e "\nContainer for ${PACKAGE} did not start yet..."
          echo -e "\nRetrying '/api/docker/containers'..."
      else
        echo "Retry $count/$retries exited, no more retries left for ${apiUrl}/docker/containers."
        exit 1
      fi
    done
    if [[ "${apiUrl}" == *"127.0.0.1"* ]]; then
      count=0
      retries=20
      until .github/scripts/check-output.sh "docker ps" "${PACKAGE}"
      do
        wait=$((2 ** count))
        count=$((count + 1))
        if [ "$count" -lt "$retries" ]; then
            sleep $wait
            echo -e "\nContainer for ${PACKAGE} did not start yet..."
            echo -e "\nRetrying 'docker ps'..."
        else
          echo "Retry $count/$retries exited, no more retries left for docker ps."
          exit 1
        fi
      done
    fi
fi

cd ${DIR} || exit 1

npm run test -- --skip-build --skip-publish --record --csv --host ${alias} --verbose

cd ${crnt} || exit 1
