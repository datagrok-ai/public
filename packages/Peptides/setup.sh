#!/bin/bash

./setup-unlink-clean.sh

GREEN='\e[0;32m'
NO_COLOR='\e[0m'

package_dir=$(pwd)

dirs=(
  "../../js-api/"
  "../../libraries/utils/"
  "../../libraries/ml/"
  "../../libraries/bio/"
)

for dir in ${dirs[@]}; do
  cd $package_dir
  cd $dir 
  echo -e $GREEN npm install in $(pwd) $NO_COLOR
  npm install 
  echo -e $GREEN npm link in $(pwd) $NO_COLOR
  npm link
done

for dir in ${dirs[@]}; do
  cd $package_dir
  cd $dir 
  if [ $dir != "../../js-api/" ]; then
    echo -e $GREEN npm link-all in $(pwd) $NO_COLOR
    npm run link-all
  fi
  echo -e $GREEN npm run build in$(pwd) $NO_COLOR
  npm run build || exit
done

cd $package_dir
npm run link-all
