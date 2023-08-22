#!/bin/bash

settings_file="runs/test_settings.conf"
num_iterations=1
selected_tests=""
default_tests=$(printf "sp3mm_test_serial\nsp3mm_test_upperbound\nsp3mm_test_upperbound_1D\nsp3mm_test_exact_rb_tree\nsp3mm_test_exact_indices\nsp3mm_test_omp\n")

# Function to display script usage/help
display_help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --with-serial            runs serial version of spmm"
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
        --with-serial)
            selected_tests=$(printf "sp3mm_test_serial\n$selected_tests")
            shift ;;
        --with-upperbound)
            selected_tests=$(printf "sp3mm_test_upperbound\n$selected_tests")
            shift ;;
        --with-upperbound-1D)
            selected_tests=$(printf "sp3mm_test_upperbound_1D\n$selected_tests")
            shift ;;
        --with-exact-rb-tree)
            selected_tests=$(printf "sp3mm_test_exact_rb_tree\n$selected_tests")
            shift ;;
        --with-exact-indices)
            selected_tests=$(printf "sp3mm_test_exact_indices\n$selected_tests")
            shift ;;
        --with-omp)
            selected_tests=$(printf "sp3mm_test_omp\n$selected_tests")
            shift ;;
        --with-omp-1d)
            selected_tests=$(printf "sp3mm_test_omp_1d\n$selected_tests")
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

if [[ $selected_tests == "" ]]; then
    echo "$default_tests" >> "$settings_file"
else
    # trimming selected tests
    echo "$selected_tests" >> "$settings_file"
fi

echo "Settings saved to $settings_file"
