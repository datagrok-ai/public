package=Chem
find packages/$package/dockerfiles -type f -name "container.json" -exec sh -c '
  for file do
    if [[ $(jq -r ".on_demand" "$file") == "true" ]]; then
      echo Skipped
    else
      container_name=$(basename $(dirname $file))
      if [[ $container_name == "dockerfiles" ]]; then
        echo $(.github/scripts/check-output.sh "docker ps" $Chem)
      else
        .github/scripts/check-output.sh "docker ps" $container_name
      fi
    fi
  done
' sh {} +
