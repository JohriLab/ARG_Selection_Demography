#!/bin/bash

#Built with ChatGPT3.5 in Sep 2023

# Check if a filename argument is provided
if [ $# -ne 1 ]; then
  echo "Usage: $0 <input_file>"
  exit 1
fi

input_file="$1"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Input file '$input_file' does not exist."
  exit 1
fi

# Loop through the lines of the input file
while read -r line; do
  # Extract the two integers from columns 1 and 2
  start=$(echo "$line" | awk '{print $1}')
  end=$(echo "$line" | awk '{print $2}')

  # Check if the input is valid (i.e., end > start)
  if [ "$end" -lt "$start" ]; then
    echo "Invalid input in line: $line. Second integer must be greater than the first."
    continue
  fi

  # Loop through and print the numbers
  for ((i = start; i <= end; i++)); do
    echo "$i"
  done
done < "$input_file"

#./print_nums.sh 2col_exons.bed > exon_positions.txt
