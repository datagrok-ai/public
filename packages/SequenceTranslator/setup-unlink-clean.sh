#!/bin/bash
package_dir=$(pwd)

GREEN='\e[0;32m'
NO_COLOR='\e[0m'

dirs=(
  "../../js-api/"
  "../../libraries/utils/"
  "../../libraries/bio/"
)

npm uninstall --location=global datagrok-api @datagrok-libraries/utils @datagrok-libraries/ml @datagrok-libraries/bio

for dir in ${dirs[@]}; do
  cd $package_dir
  cd $dir 
  echo -e $GREEN Removing node_modules and dist in $(pwd) $NO_COLOR
  rm -rf node_modules dist
  # rm package-lock.json
done
