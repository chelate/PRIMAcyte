#!/bin/bash

# Get the directory where the script is located
SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"

# Define the full path to the data directory
DATA_DIR="$SCRIPT_DIR/data"

# Define the absolute path to the home directory
HOME_DIR="$(echo ~)"

# Check if the data directory exists, create if not
if [ ! -d "$DATA_DIR" ]; then
    mkdir -p "$DATA_DIR"
fi

# Define the Docker run command
DOCKER_CMD="docker run -v $DATA_DIR:/data \
                      --platform linux/amd64 \
                      -u $(id -u):$(id -g) \
                      -p 8888:8888 \
                      -v /tmp/.X11-unix:/tmp/.X11-unix \
                      -v $HOME_DIR/.Xauthority:/home/steinbock/.Xauthority:ro \
                      -e DISPLAY=host.docker.internal:0 \
                      ghcr.io/bodenmillergroup/steinbock:0.16.1-cellpose"

# Set the alias for the current shell session
alias steinbock="$DOCKER_CMD"

echo "Alias 'steinbock' set to docker steinbock:0.16.1-cellpose" 
echo "mount $DATA_DIR as /data"
#$DOCKER_CMD imc unzip

#echo "Executing 'steinbock panel'..."
#$DOCKER_CMD panel

pluto() {
  PLUTO_ARGS="dismiss_update_notification=true, require_secret_for_open_links=true, require_secret_for_access=true"
  case "$2" in
      yes) REVISE="using Revise;" ;;
      *)   REVISE="" ;;
  esac
  if [[ -d "$HOME/Projects/Pluto.jl" ]]; then
    PLUTO_PROJECT=$HOME/Projects/Pluto.jl
  else
    PLUTO_PROJECT=@pluto
  fi
  julia --threads=auto --project=$PLUTO_PROJECT -e "
    $REVISE
    import Pluto;
    notebook_to_open = \"$1\"
    if length(notebook_to_open) == 0
      Pluto.run($PLUTO_ARGS)
    else
      @info \"opening notebook\" notebook_to_open
      Pluto.run(notebook=notebook_to_open, $PLUTO_ARGS)
    end
  "
}