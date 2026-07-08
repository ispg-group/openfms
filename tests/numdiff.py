#!/usr/bin/env python3
# /// script
# requires-python = ">=3.6"
# dependencies = []
# ///
import argparse
import decimal
import os.path
import sys

# This scripts is compares differences between files
# while disregarding insignificant numerical differences.

# This script parses output from command:
# diff -y file1 file2

# It returns exit code 1 when difference is found and 0 otherwise

decimal.getcontext().prec = 17
# decimal.getcontext().rounding = "ROUND_DOWN"


def read_cmd():
    """Function for reading command line options.
    returns the name of the input file."""
    desc = """Compare numerical differences between test and reference files.
            Parses the output of 'diff -y file file.ref'"""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("input_file", metavar="INPUT_FILE", help="input file")
    opts = parser.parse_args()
    return opts.input_file


def read_custom_threshold(fname):
    """Reads a custom absolute threshold exponent from a file in a test folder"""
    with open(fname, encoding="utf-8") as f:
        thr = int(f.read().strip())

    if thr < 2 or thr > 20:
        print(f"ERROR: Threshold 10^-{thr} read from {fname} is out of range")
        sys.exit(1)

    return 10 ** (-thr)


def compare_numbers(test, reference, absolute_tolerance):
    delta = abs(float(reference) - float(test))
    if delta > absolute_tolerance:
        print("Significant numerical difference!")
        print(f"Reference = {reference} Test = {test} Delta = {delta:.2e}")
        sys.exit(1)


def compare_lines(numbers1, numbers2, absolute_tolerance):
    """Compare each element of numbers1 and numbers2"""
    for test, reference in zip(numbers1, numbers2):
        if test == reference:
            continue

        try:
            # Rounding to significant precision
            dec1 = +decimal.Decimal(test)
            dec2 = +decimal.Decimal(reference)
        except decimal.InvalidOperation:
            # workaround for Last_Bundle.txt, where we have a string like:
            # (0.98831589475674242,-1.79586289170714135E-003)
            new_test = [s.strip('()') for s in test.split(',')]
            new_reference = [s.strip('()') for s in reference.split(',')]
            if len(new_test) == 2 and len(new_reference) == 2:
                compare_lines(new_test, new_reference, absolute_tolerance)
                continue
            else:
                print("Non-numerical difference found!")
                print(f"Test     :{test}")
                print(f"Reference:{reference}")
                sys.exit(1)

        if dec1 == dec2:
            return

        # Compare sign, exponents, and digits separately
        sign1, digits1, exp1 = dec1.as_tuple()
        sign2, digits2, exp2 = dec2.as_tuple()
        if sign1 != sign2 or exp1 != exp2:
            compare_numbers(test, reference, absolute_tolerance)
            continue

        ints1 = ""
        ints2 = ""
        # Compare last 10 significant digits
        for j in digits1[-10:]:
            ints1 += str(j)
        for j in digits2[-10:]:
            ints2 += str(j)
        ints1 = int(ints1)
        ints2 = int(ints2)
        # This allows to get rid off insignificant rounding errors
        if abs(ints1 - ints2) != 1:
            compare_numbers(test, reference, absolute_tolerance)


def parse_diff(fname, absolute_tolerance):
    with open(fname, encoding="locale") as f:
        # If the file is empty something is wrong
        # (e.g. file for comparison was not even generated)
        if len(f.read()) == 0:
            print(f"File '{fname}' is empty!")
            sys.exit(1)
        f.seek(0)

        for line in f:
            split = line.split()
            # Exit if there are non-numerical differences
            if "<" in split or ">" in split:
                print("There's an extra line!")
                print(line)
                sys.exit(1)
            # Skip line if there is no difference
            if "|" not in split:
                continue

            split = line.split("|")
            diff1 = split[0].split()
            diff2 = split[1].split()
            if len(diff1) != len(diff2):
                print("Number of columns differ!")
                print("Expected:")
                print(split[1])
                print("Test:")
                print(split[0])
                sys.exit(1)
            compare_lines(diff1, diff2, absolute_tolerance)


if __name__ == "__main__":
    absolute_tolerance = 2e-15
    # If an individual test needs a higher tolerance, it can do so
    # by specifying the exponent in file "NUM_THRE" in the test folder.
    # For example, NUM_THRE containing "3" will set tolerance = 1e-3
    THR_FNAME = "NUM_THRE"
    if os.path.isfile(THR_FNAME):
        absolute_tolerance = read_custom_threshold(THR_FNAME)

    inpfile = read_cmd()
    print("\nComparing numerical differences in file " + inpfile)
    parse_diff(inpfile, absolute_tolerance)
    sys.exit(0)
