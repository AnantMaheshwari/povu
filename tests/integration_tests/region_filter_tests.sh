#!/bin/bash
#
# Integration tests for genomic region filtering feature
# Tests the --region flag for the call subcommand
#

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
POVU_BIN="$REPO_ROOT/bin/povu"
TEST_DATA_DIR="$REPO_ROOT/tests/data"
TEST_OUTPUT_DIR="/tmp/povu_region_tests_$$"

# Test status counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Create test output directory
mkdir -p "$TEST_OUTPUT_DIR"

# Cleanup function
cleanup() {
    rm -rf "$TEST_OUTPUT_DIR"
}
trap cleanup EXIT

# Helper functions
print_test_header() {
    echo ""
    echo "======================================================================"
    echo "TEST: $1"
    echo "======================================================================"
}

print_success() {
    echo -e "${GREEN}✓ PASS${NC}: $1"
    ((TESTS_PASSED++))
}

print_failure() {
    echo -e "${RED}✗ FAIL${NC}: $1"
    ((TESTS_FAILED++))
}

print_info() {
    echo -e "${YELLOW}INFO${NC}: $1"
}

run_test() {
    ((TESTS_RUN++))
}

# ============================================================================
# Test 1: Verify --region flag appears in help
# ============================================================================
test_region_flag_in_help() {
    print_test_header "Verify --region flag appears in help"
    run_test
    
    if timeout 5 "$POVU_BIN" call --help 2>&1 | grep -q "\-\-region"; then
        print_success "Region flag appears in help text"
    else
        print_failure "Region flag missing from help text"
    fi
}

# ============================================================================
# Test 2: Test invalid region format handling
# ============================================================================
test_invalid_region_format() {
    print_test_header "Test invalid region format handling"
    run_test
    
    # Test with a simple synthetic GFA (no decomposition needed)
    local TEST_GFA="$TEST_OUTPUT_DIR/simple.gfa"
    
    # Create a minimal GFA file
    cat > "$TEST_GFA" << 'EOF'
H	VN:Z:1.0
S	1	ACGT
S	2	TTGG
L	1	+	2	+	0M
W	sample1	0	ref1	0	4	>1>2
EOF
    
    print_info "Testing region with missing colon..."
    if timeout 10 "$POVU_BIN" call -i "$TEST_GFA" -P ref1 --region "invalidregion" --stdout 2>&1 | grep -q "Invalid region format.*missing colon"; then
        print_success "Correctly detected missing colon in region format"
    else
        print_failure "Did not detect missing colon in region format"
    fi
    
    run_test
    print_info "Testing region with missing dash..."
    if timeout 10 "$POVU_BIN" call -i "$TEST_GFA" -P ref1 --region "ref1:100" --stdout 2>&1 | grep -q "Invalid region format.*missing dash"; then
        print_success "Correctly detected missing dash in region format"
    else
        print_failure "Did not detect missing dash in region format"
    fi
    
    run_test
    print_info "Testing region with non-numeric positions..."
    if timeout 10 "$POVU_BIN" call -i "$TEST_GFA" -P ref1 --region "ref1:abc-def" --stdout 2>&1 | grep -q "Invalid region format.*non-numeric"; then
        print_success "Correctly detected non-numeric positions"
    else
        print_failure "Did not detect non-numeric positions"
    fi
    
    run_test
    print_info "Testing region with start >= end..."
    if timeout 10 "$POVU_BIN" call -i "$TEST_GFA" -P ref1 --region "ref1:100-50" --stdout 2>&1 | grep -q "start position.*must be less than end"; then
        print_success "Correctly detected start >= end"
    else
        print_failure "Did not detect start >= end"
    fi
}

# ============================================================================
# Test 3: Test backward compatibility (no region specified)
# ============================================================================
test_backward_compatibility() {
    print_test_header "Test backward compatibility (no region specified)"
    run_test
    
    local TEST_GFA="$TEST_OUTPUT_DIR/simple2.gfa"
    
    # Create a minimal GFA file with reference
    cat > "$TEST_GFA" << 'EOF'
H	VN:Z:1.0
S	1	ACGT
S	2	TTGG
S	3	AAAA
L	1	+	2	+	0M
L	2	+	3	+	0M
W	sample1	0	chr1	0	12	>1>2>3
EOF
    
    print_info "Running decompose to generate PVST..."
    "$POVU_BIN" decompose -i "$TEST_GFA" -o "$TEST_OUTPUT_DIR" >/dev/null 2>&1 || true
    
    print_info "Running call without region (backward compatibility)..."
    if "$POVU_BIN" call -i "$TEST_GFA" -f "$TEST_OUTPUT_DIR" -P chr1 --stdout 2>&1 | grep -q "^##fileformat=VCFv4.2"; then
        print_success "Call works without region parameter (backward compatible)"
    else
        print_failure "Call failed without region parameter"
    fi
}

# ============================================================================
# Test 4: Test with valid region on synthetic data
# ============================================================================
test_valid_region_synthetic() {
    print_test_header "Test with valid region on synthetic data"
    run_test
    
    local TEST_GFA="$TEST_OUTPUT_DIR/variant_graph.gfa"
    
    # Create a GFA with some variation
    cat > "$TEST_GFA" << 'EOF'
H	VN:Z:1.0
S	1	AAAA
S	2	TTTT
S	3	GGGG
S	4	CCCC
S	5	AAAA
S	6	TTTT
L	1	+	2	+	0M
L	1	+	3	+	0M
L	2	+	4	+	0M
L	3	+	4	+	0M
L	4	+	5	+	0M
L	4	+	6	+	0M
L	5	+	7	+	0M
L	6	+	7	+	0M
S	7	GGGG
W	sample1	0	chr1	0	24	>1>2>4>5>7
W	sample2	0	chr1	0	24	>1>3>4>6>7
EOF
    
    print_info "Running decompose to generate PVST..."
    "$POVU_BIN" decompose -i "$TEST_GFA" -o "$TEST_OUTPUT_DIR" >/dev/null 2>&1 || true
    
    local VCF_NO_REGION="$TEST_OUTPUT_DIR/no_region.vcf"
    local VCF_WITH_REGION="$TEST_OUTPUT_DIR/with_region.vcf"
    
    print_info "Running call without region..."
    "$POVU_BIN" call -i "$TEST_GFA" -f "$TEST_OUTPUT_DIR" -P chr1 --stdout > "$VCF_NO_REGION" 2>&1 || true
    
    print_info "Running call with region chr1:0-12..."
    "$POVU_BIN" call -i "$TEST_GFA" -f "$TEST_OUTPUT_DIR" -P chr1 --region "chr1:0-12" --stdout > "$VCF_WITH_REGION" 2>&1 || true
    
    # Check that both generated VCF headers
    if grep -q "^##fileformat=VCFv4.2" "$VCF_NO_REGION"; then
        print_success "Generated VCF without region filter"
    else
        print_failure "Failed to generate VCF without region filter"
    fi
    
    run_test
    if grep -q "^##fileformat=VCFv4.2" "$VCF_WITH_REGION"; then
        print_success "Generated VCF with region filter"
    else
        print_failure "Failed to generate VCF with region filter"
    fi
    
    # Count variants
    run_test
    local COUNT_NO_REGION=$(grep -c "^chr1" "$VCF_NO_REGION" || echo 0)
    local COUNT_WITH_REGION=$(grep -c "^chr1" "$VCF_WITH_REGION" || echo 0)
    
    print_info "Variants without region: $COUNT_NO_REGION"
    print_info "Variants with region chr1:0-12: $COUNT_WITH_REGION"
    
    if [ "$COUNT_WITH_REGION" -le "$COUNT_NO_REGION" ]; then
        print_success "Region filtering reduced or maintained variant count (expected behavior)"
    else
        print_failure "Region filtering increased variant count (unexpected)"
    fi
}

# ============================================================================
# Test 5: Test with region outside graph bounds
# ============================================================================
test_region_outside_bounds() {
    print_test_header "Test with region outside graph bounds"
    run_test
    
    local TEST_GFA="$TEST_OUTPUT_DIR/bounds_test.gfa"
    
    # Create a minimal GFA
    cat > "$TEST_GFA" << 'EOF'
H	VN:Z:1.0
S	1	ACGT
S	2	TTGG
L	1	+	2	+	0M
W	sample1	0	chr1	0	8	>1>2
EOF
    
    print_info "Running decompose..."
    "$POVU_BIN" decompose -i "$TEST_GFA" -o "$TEST_OUTPUT_DIR" >/dev/null 2>&1 || true
    
    print_info "Testing region completely outside bounds (chr1:1000-2000)..."
    local VCF_OUT="$TEST_OUTPUT_DIR/outside_bounds.vcf"
    "$POVU_BIN" call -i "$TEST_GFA" -f "$TEST_OUTPUT_DIR" -P chr1 --region "chr1:1000-2000" --stdout > "$VCF_OUT" 2>&1 || true
    
    # Should generate empty or minimal VCF (header only)
    if grep -q "^##fileformat=VCFv4.2" "$VCF_OUT"; then
        local VAR_COUNT=$(grep -c "^chr1" "$VCF_OUT" || echo 0)
        if [ "$VAR_COUNT" -eq 0 ]; then
            print_success "Region outside bounds correctly returns no variants"
        else
            print_failure "Region outside bounds returned $VAR_COUNT variants (expected 0)"
        fi
    else
        print_failure "Failed to generate VCF for region outside bounds"
    fi
}

# ============================================================================
# Test 6: Test with non-existent reference path
# ============================================================================
test_nonexistent_reference() {
    print_test_header "Test with non-existent reference path"
    run_test
    
    local TEST_GFA="$TEST_OUTPUT_DIR/ref_test.gfa"
    
    cat > "$TEST_GFA" << 'EOF'
H	VN:Z:1.0
S	1	ACGT
S	2	TTGG
L	1	+	2	+	0M
W	sample1	0	chr1	0	8	>1>2
EOF
    
    print_info "Running decompose..."
    "$POVU_BIN" decompose -i "$TEST_GFA" -o "$TEST_OUTPUT_DIR" >/dev/null 2>&1 || true
    
    print_info "Testing with non-existent reference 'chr99'..."
    if "$POVU_BIN" call -i "$TEST_GFA" -f "$TEST_OUTPUT_DIR" -P chr1 --region "chr99:0-100" --stdout 2>&1 | grep -q "Reference path.*not found"; then
        print_success "Correctly detected non-existent reference path"
    else
        print_failure "Did not detect non-existent reference path"
    fi
}

# ============================================================================
# Main test execution
# ============================================================================
main() {
    echo ""
    echo "======================================================================"
    echo "Genomic Region Filtering Integration Tests"
    echo "======================================================================"
    echo ""
    
    # Check if povu binary exists
    if [ ! -f "$POVU_BIN" ]; then
        echo -e "${RED}ERROR${NC}: povu binary not found at $POVU_BIN"
        echo "Please build povu first: cmake -H. -Bbuild && cmake --build build -- -j 3"
        exit 1
    fi
    
    print_info "Using povu binary: $POVU_BIN"
    print_info "Test data directory: $TEST_DATA_DIR"
    print_info "Test output directory: $TEST_OUTPUT_DIR"
    echo ""
    
    # Run all tests
    test_region_flag_in_help
    test_invalid_region_format
    test_backward_compatibility
    test_valid_region_synthetic
    test_region_outside_bounds
    test_nonexistent_reference
    
    # Print summary
    echo ""
    echo "======================================================================"
    echo "Test Summary"
    echo "======================================================================"
    echo "Tests run:    $TESTS_RUN"
    echo -e "Tests passed: ${GREEN}$TESTS_PASSED${NC}"
    if [ $TESTS_FAILED -gt 0 ]; then
        echo -e "Tests failed: ${RED}$TESTS_FAILED${NC}"
        exit 1
    else
        echo -e "Tests failed: ${GREEN}0${NC}"
        echo ""
        echo -e "${GREEN}All tests passed!${NC}"
        exit 0
    fi
}

main "$@"
