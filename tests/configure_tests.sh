#!/bin/bash

settings_file="runs/test_settings.conf"
num_iterations=1
nb_selected_tests=0
selected_tests=""
nb_default_tests=12
default_tests="serial \
            serial_reversed \
            omp_gustavson \
            omp_gustavson_reversed \
            omp_gustavson_1d\
            omp_gustavson_1d_reversed \
            serial_rb_tree \
            serial_rb_tree_reversed \
            omp_rb_tree \
            omp_rb_tree_reversed \
            omp_two_pass \
            omp_two_pass_reversed"

# Function to display script usage/help
display_help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --with-serial            runs serial version of spmm"
    echo "  --with-serial_reversed   runs serial version of spmm but order of the multiplications is reversed"
    echo "  --with-upperbound        runs basic upperbounded version"
    echo "  --with-upperbound-1D     runs 1D version of the upperbounded algorithm"
    echo "  --with-exact-rb-tree     runs exact memory allocation version with red-black-tree distribution"
    echo "  --with-exact-indices     runs exact memory allocation version with indices computed first"
    echo "  -i, --iterations <num>   Set number of iterations (default: 1)"
    echo "  -h, --help               Display this help"
}

# Parse command-line options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --with-*)
            nb_selected_tests=$(( $nb_selected_tests + 1 ))
            selected_tests="$selected_tests $(echo $1 | cut -d '-' -f 4-)"
            shift ;;
        -i | --iterations)
            if [[ -n "$2" && "$2" =~ ^[0-9]+$ ]]; then
                num_iterations="$2"
                shift 2
            else
                echo "Invalid argument for -i/--iterations. Expected a number."
                exit 1
            fi ;;
        -h | --help)
            display_help
            exit 0 ;;
        *)
            echo "Unknown option: $1"
            display_help
            exit 1 ;;
    esac
done

# Add settings
echo "$num_iterations" > "$settings_file"
echo ".true." >> "$settings_file"

if [[ $selected_tests == "" ]]; then
    echo $nb_default_tests >> "$settings_file"
    for test in $default_tests
    do
        echo $test >> "$settings_file"
    done
else
    echo $nb_selected_tests >> "$settings_file"
    for test in $selected_tests
    do
        echo $test >> "$settings_file"
    done
fi

echo "Settings saved to $settings_file"
