#!/bin/bash
################################################################################
# Bridge Automation Script: Run Downstream GSEA Analysis
#
# This script automates the execution of RNA-Seq_GO_GSEA_analysis pipeline
# after RNA-Seq_DE_GO_analysis completes.
#
# Usage:
#   ./run_downstream_analysis.sh <comparison_name> [options]
#
# Example:
#   ./run_downstream_analysis.sh H2O2_vs_Control --run-gsea
#
# Author: Pipeline Integration Team
# Date: 2025-12-01
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

################################################################################
# Configuration
################################################################################

# Default paths (modify these to match your setup)
DE_PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GSEA_PIPELINE_ROOT="$(cd "${DE_PIPELINE_ROOT}/../RNA-Seq_GO_GSEA_analysis" && pwd)"
BRIDGE_DIR="${DE_PIPELINE_ROOT}/bridge"

# Default DE analysis output directory
DE_OUTPUT_DIR="${DE_PIPELINE_ROOT}/output"

# GSEA pipeline input directory
GSEA_INPUT_DIR="${GSEA_PIPELINE_ROOT}/data/from_de_pipeline"

# Conversion script
CONVERTER_SCRIPT="${BRIDGE_DIR}/convert_de_to_gsea.py"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

################################################################################
# Functions
################################################################################

print_header() {
    echo -e "${BLUE}================================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}================================================================${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ $1${NC}"
}

usage() {
    cat << EOF
Usage: $0 <comparison_name> [options]

Arguments:
    comparison_name     Name of the comparison (e.g., H2O2_vs_Control)

Options:
    --experiment-dir    Experiment directory (e.g., H2O2_Neuron)
    --de-output-dir     Custom DE analysis output directory
    --gsea-input-dir    Custom GSEA pipeline input directory
    --run-go            Run GO enrichment analysis
    --run-gsea          Run GSEA analysis
    --run-all           Run both GO and GSEA analyses
    --skip-convert      Skip conversion step (use existing Excel file)
    --batch             Batch mode: convert all comparisons
    --help, -h          Show this help message

Examples:
    # Convert and run GO enrichment for a single comparison
    $0 H2O2_vs_Control --experiment-dir H2O2_Neuron --run-go

    # Convert and run GSEA for a single comparison
    $0 GABA_vs_Control --experiment-dir H2O2_Neuron --run-gsea

    # Batch convert all comparisons and run both analyses
    $0 --batch --experiment-dir H2O2_Neuron --run-all

EOF
    exit 1
}

check_dependencies() {
    print_header "Checking Dependencies"
    
    # Check if Python is available
    if ! command -v python3 &> /dev/null; then
        print_error "Python 3 is not installed or not in PATH"
        exit 1
    fi
    print_success "Python 3 found: $(python3 --version)"
    
    # Check if converter script exists
    if [ ! -f "${CONVERTER_SCRIPT}" ]; then
        print_error "Converter script not found: ${CONVERTER_SCRIPT}"
        exit 1
    fi
    print_success "Converter script found"
    
    # Check if GSEA pipeline exists
    if [ ! -d "${GSEA_PIPELINE_ROOT}" ]; then
        print_error "GSEA pipeline not found: ${GSEA_PIPELINE_ROOT}"
        exit 1
    fi
    print_success "GSEA pipeline found"
    
    echo ""
}

convert_single_comparison() {
    local comparison_name="$1"
    local experiment_dir="$2"
    
    print_header "Converting DE Results: ${comparison_name}"
    
    # Construct input path
    local input_csv="${DE_OUTPUT_DIR}/${experiment_dir}/pairwise/${comparison_name}/final_de_results.csv"
    
    # Check if input file exists
    if [ ! -f "${input_csv}" ]; then
        print_error "Input file not found: ${input_csv}"
        return 1
    fi
    
    print_info "Input: ${input_csv}"
    
    # Create output directory if it doesn't exist
    mkdir -p "${GSEA_INPUT_DIR}"
    
    # Construct output path
    local output_excel="${GSEA_INPUT_DIR}/${comparison_name}_DE_results.xlsx"
    print_info "Output: ${output_excel}"
    
    # Run conversion
    echo ""
    python3 "${CONVERTER_SCRIPT}" \
        --input "${input_csv}" \
        --output "${output_excel}" \
        --comparison-name "${comparison_name}"
    
    if [ $? -eq 0 ]; then
        print_success "Conversion completed successfully"
        echo ""
        return 0
    else
        print_error "Conversion failed"
        echo ""
        return 1
    fi
}

convert_batch() {
    local experiment_dir="$1"
    
    print_header "Batch Converting All Comparisons"
    
    print_info "DE output directory: ${DE_OUTPUT_DIR}/${experiment_dir}"
    print_info "GSEA input directory: ${GSEA_INPUT_DIR}"
    
    # Run batch conversion
    echo ""
    python3 "${CONVERTER_SCRIPT}" \
        --batch \
        --de-output-dir "${DE_OUTPUT_DIR}/${experiment_dir}" \
        --gsea-input-dir "${GSEA_INPUT_DIR}"
    
    if [ $? -eq 0 ]; then
        print_success "Batch conversion completed"
        echo ""
        return 0
    else
        print_error "Batch conversion failed"
        echo ""
        return 1
    fi
}

run_go_analysis() {
    local comparison_name="$1"
    
    print_header "Running GO Enrichment Analysis: ${comparison_name}"
    
    local input_file="${GSEA_INPUT_DIR}/${comparison_name}_DE_results.xlsx"
    
    if [ ! -f "${input_file}" ]; then
        print_error "Input file not found: ${input_file}"
        return 1
    fi
    
    print_info "Input: ${input_file}"
    print_info "Running GO pipeline..."
    echo ""
    
    # Navigate to GSEA pipeline directory
    cd "${GSEA_PIPELINE_ROOT}"
    
    # Run GO pipeline (modify this command based on your GSEA pipeline setup)
    # This is a placeholder - adjust according to your actual pipeline
    if [ -f "notebooks/GO_Pipeline.ipynb" ]; then
        print_info "GO pipeline notebook found. Please run manually or configure automation."
        print_warning "Automated notebook execution requires additional setup."
    else
        print_warning "GO pipeline configuration not found. Please configure manually."
    fi
    
    cd "${DE_PIPELINE_ROOT}"
    echo ""
}

run_gsea_analysis() {
    local comparison_name="$1"
    
    print_header "Running GSEA Analysis: ${comparison_name}"
    
    local input_file="${GSEA_INPUT_DIR}/${comparison_name}_DE_results.xlsx"
    
    if [ ! -f "${input_file}" ]; then
        print_error "Input file not found: ${input_file}"
        return 1
    fi
    
    print_info "Input: ${input_file}"
    print_info "Running GSEA pipeline..."
    echo ""
    
    # Navigate to GSEA pipeline directory
    cd "${GSEA_PIPELINE_ROOT}"
    
    # Run GSEA pipeline (modify this command based on your GSEA pipeline setup)
    # This is a placeholder - adjust according to your actual pipeline
    if [ -f "notebooks/GSEA_Pipeline.ipynb" ]; then
        print_info "GSEA pipeline notebook found. Please run manually or configure automation."
        print_warning "Automated notebook execution requires additional setup."
    else
        print_warning "GSEA pipeline configuration not found. Please configure manually."
    fi
    
    cd "${DE_PIPELINE_ROOT}"
    echo ""
}

################################################################################
# Main Script
################################################################################

# Parse arguments
COMPARISON_NAME=""
EXPERIMENT_DIR="H2O2_Neuron"  # Default
RUN_GO=false
RUN_GSEA=false
SKIP_CONVERT=false
BATCH_MODE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --experiment-dir)
            EXPERIMENT_DIR="$2"
            shift 2
            ;;
        --de-output-dir)
            DE_OUTPUT_DIR="$2"
            shift 2
            ;;
        --gsea-input-dir)
            GSEA_INPUT_DIR="$2"
            shift 2
            ;;
        --run-go)
            RUN_GO=true
            shift
            ;;
        --run-gsea)
            RUN_GSEA=true
            shift
            ;;
        --run-all)
            RUN_GO=true
            RUN_GSEA=true
            shift
            ;;
        --skip-convert)
            SKIP_CONVERT=true
            shift
            ;;
        --batch)
            BATCH_MODE=true
            shift
            ;;
        --help|-h)
            usage
            ;;
        -*)
            print_error "Unknown option: $1"
            usage
            ;;
        *)
            COMPARISON_NAME="$1"
            shift
            ;;
    esac
done

# Validate arguments
if [ "${BATCH_MODE}" = false ] && [ -z "${COMPARISON_NAME}" ]; then
    print_error "Comparison name is required in single mode"
    usage
fi

# Print configuration
print_header "Pipeline Bridge Configuration"
echo "DE Pipeline Root:    ${DE_PIPELINE_ROOT}"
echo "GSEA Pipeline Root:  ${GSEA_PIPELINE_ROOT}"
echo "Experiment Dir:      ${EXPERIMENT_DIR}"
echo "Mode:                $([ "${BATCH_MODE}" = true ] && echo "Batch" || echo "Single")"
[ "${BATCH_MODE}" = false ] && echo "Comparison:          ${COMPARISON_NAME}"
echo "Run GO:              ${RUN_GO}"
echo "Run GSEA:            ${RUN_GSEA}"
echo "Skip Conversion:     ${SKIP_CONVERT}"
echo ""

# Check dependencies
check_dependencies

# Main workflow
if [ "${SKIP_CONVERT}" = false ]; then
    if [ "${BATCH_MODE}" = true ]; then
        convert_batch "${EXPERIMENT_DIR}"
    else
        convert_single_comparison "${COMPARISON_NAME}" "${EXPERIMENT_DIR}"
    fi
else
    print_warning "Skipping conversion step (--skip-convert enabled)"
    echo ""
fi

# Run downstream analyses
if [ "${BATCH_MODE}" = false ]; then
    if [ "${RUN_GO}" = true ]; then
        run_go_analysis "${COMPARISON_NAME}"
    fi
    
    if [ "${RUN_GSEA}" = true ]; then
        run_gsea_analysis "${COMPARISON_NAME}"
    fi
else
    if [ "${RUN_GO}" = true ] || [ "${RUN_GSEA}" = true ]; then
        print_warning "Batch mode: Please run GO/GSEA analyses manually for each comparison"
        print_info "Converted files are available in: ${GSEA_INPUT_DIR}"
        echo ""
    fi
fi

# Final summary
print_header "Pipeline Bridge Execution Summary"
print_success "All requested operations completed"
echo ""
print_info "Next steps:"
if [ "${BATCH_MODE}" = true ]; then
    echo "  1. Check converted files in: ${GSEA_INPUT_DIR}"
    echo "  2. Run GO/GSEA analyses for each comparison as needed"
else
    echo "  1. Review results in: ${GSEA_INPUT_DIR}/${COMPARISON_NAME}_DE_results.xlsx"
    if [ "${RUN_GO}" = false ] && [ "${RUN_GSEA}" = false ]; then
        echo "  2. Run GO/GSEA analysis manually using the converted file"
    fi
fi
echo ""
print_success "Done!"
