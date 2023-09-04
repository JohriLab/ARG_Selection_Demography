import sys

if len(sys.argv) != 3:
    print("Usage: python script.py file1.txt file2.txt")
    sys.exit(1)

file1_path = sys.argv[1]
file2_path = sys.argv[2]

# Read the content of file 1 and extract the two columns as ranges
with open(file1_path, 'r') as file1:
    lines = file1.readlines()
    ranges = [(int(line.split()[1]), int(line.split()[2])) for line in lines]

# Read the content of file 2 and extract the integers
with open(file2_path, 'r') as file2:
    integers = [int(line.strip()) for line in file2.readlines()]

# Iterate through the integers and check if they are within any of the ranges
for num in integers:
    for range_start, range_end in ranges:
        if range_start <= num < range_end:
            print(num)
            break  # No need to check other ranges if the number is already found within one

