#!/bin/bash

# ANSI escape codes for colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# Check if username is provided as an argument
if [ $# -eq 0 ]; then
  echo "Usage: $0 <username>"
  exit 1
fi

# Define the username
username=$1

# Get the current host's name
current_host=$(hostname)

# Define the list of hosts with the local host moved to the top
hosts=(
  "$current_host"  # local host
  "MI104-01"
  "MI104-02"
  "MI104-03"
  "MI104-04"
  "MI104-05"
  "MI104-06"
  "MI104-07"
  "MI104-08"
  "MI104-09"
  "MI104-10"
  "MI104-11"
  "MI104-12"
  "MI104-13"
  "MI104-14"
  "MI104-15"
  "MI104-16"
  "MI104-17"
  "MI104-18"
  "MI104-19"
  "MI104-20"
  "MI105-14"
  "MI105-13"
  "MI105-12"
  "MI105-11"
  "MI105-10"
  "MI105-09"
  "MI105-08"
  "MI105-07"
  "MI105-06"
  "MI105-05"
  "MI105-04"
  "MI105-03"
  "MI105-02"
  "MI105-01"
  "MI121-01"
  "MI121-02"
  "MI121-03"
  "MI121-04"
  "MI121-05"
  "MI121-06"
  "MI121-07"
  "MI121-08"
)

# Remove the existing "hosts" file
rm -f hosts

# Loop through each host
for host in "${hosts[@]}"; do
  # Skip the local host during testing if it's not the first iteration
  if [ "$host" == "$current_host" ] && [ "$tested_local_host" == "true" ]; then
    continue
  fi

  echo -n "Trying $host... "

  # Connect to the host using SSH and retrieve the number of threads
  threads=$(ssh -o ConnectTimeout=2 -o BatchMode=yes "$username@$host" nproc --all 2>/dev/null)

  # Check if the SSH connection was successful and the number of threads is available
  if [ $? -eq 0 ] && [ -n "$threads" ]; then
    # Print the host information to the "hosts" file in green
    echo -e "${GREEN}Success${NC}"
    echo -e "$host slots=$threads" >> hosts

    # If the tested host is the local host, mark it as tested
    [ "$host" == "$current_host" ] && tested_local_host="true"
  else
    # Print the host information to the "hosts" file in red
    echo -e "${RED}Failed${NC}"
  fi
done


