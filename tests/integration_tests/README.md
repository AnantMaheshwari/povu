# Integration Tests for Genomic Region Filtering

This directory contains integration tests for the `--region` flag added to the `call` subcommand.

## Test Coverage

The integration tests validate:

1. **CLI Integration**: Region flag appears in help text
2. **Input Validation**: 
   - Invalid region formats (missing colon, missing dash)
   - Non-numeric positions
   - Start position >= end position
   - Empty position strings
3. **Reference Validation**: Non-existent reference paths
4. **Region Filtering**: Correct filtering of RoVs based on genomic coordinates
5. **Backward Compatibility**: Call works without region parameter
6. **Edge Cases**: Regions outside graph bounds

## Running the Tests

### Quick Test - Verify Help Text

```bash
./bin/povu call --help | grep region
```

Expected output:
```
-g[region], --region=[region]     Restrict variant calling to a genomic
                                  region (format: ref:start-end)
```

### Full Integration Test Suite

```bash
bash tests/integration_tests/region_filter_tests.sh
```

This runs all integration tests and reports pass/fail status.

## Manual Testing Examples

### Example 1: Invalid Region Format

```bash
# Create test GFA
cat > test.gfa << 'EOF'
H	VN:Z:1.1
S	1	ACGT
S	2	TTGG
L	1	+	2	+	0M
W	sample1	0	chr1	0	8	>1>2
EOF

# Test with invalid format (missing colon)
./bin/povu decompose -i test.gfa -o .
./bin/povu call -i test.gfa -f . -P chr1 --region "invalidformat" --stdout

# Expected: Error message about invalid format
```

### Example 2: Valid Region Filter

```bash
# Create test GFA with variants
cat > variant.gfa << 'EOF'
H	VN:Z:1.1
S	1	AAAA
S	2	TTTT
S	3	GGGG
S	4	CCCC
L	1	+	2	+	0M
L	1	+	3	+	0M
L	2	+	4	+	0M
L	3	+	4	+	0M
W	sample1	0	chr1	0	16	>1>2>4
W	sample2	0	chr1	0	16	>1>3>4
EOF

# Decompose graph
./bin/povu decompose -i variant.gfa -o .

# Call all variants (no region)
./bin/povu call -i variant.gfa -f . -P chr1 --stdout > all_variants.vcf

# Call variants in region chr1:0-8
./bin/povu call -i variant.gfa -f . -P chr1 --region "chr1:0-8" --stdout > region_variants.vcf

# Compare variant counts
echo "All variants: $(grep -c "^chr1" all_variants.vcf || echo 0)"
echo "Region variants: $(grep -c "^chr1" region_variants.vcf || echo 0)"
```

### Example 3: Region Outside Bounds

```bash
# Using same GFA from Example 2

# Call with region far outside graph
./bin/povu call -i variant.gfa -f . -P chr1 --region "chr1:10000-20000" --stdout > outside.vcf

# Should return header but no variants
grep -c "^chr1" outside.vcf || echo "0 variants (expected)"
```

## Test Data

The integration tests use synthetic GFA files created on-the-fly to ensure:
- Known graph structure
- Predictable variant positions
- Control over reference paths
- No external dependencies

## Expected Behavior

### Valid Region
- Filters PVST vertices based on position overlap
- Returns only RoVs overlapping the specified region
- Maintains VCF format and headers
- Logs filtering information

### Invalid Region
- Validates format (ref:start-end)
- Checks for non-numeric positions
- Ensures start < end
- Provides clear error messages

### Backward Compatibility
- Works identically without --region flag
- No changes to existing behavior
- Optional parameter

## Troubleshooting

If tests fail:

1. **Build Issues**: Ensure povu is built with latest changes
   ```bash
   cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
   cmake --build build -- -j 3
   ```

2. **Missing PVST Files**: The `call` command requires decomposed graph
   ```bash
   ./bin/povu decompose -i input.gfa -o output_dir
   ./bin/povu call -i input.gfa -f output_dir -P ref --region "ref:start-end"
   ```

3. **GFA Format**: Ensure GFA uses version 1.1 with W-lines for references
   ```
   H	VN:Z:1.1
   W	sample	haplotype	ref_name	start	end	walk
   ```

## Implementation Notes

The region filtering:
- Uses step indices from `vertex_to_step_matrix_` for position lookups
- Filters at PVST vertex level before walk generation
- Checks source and sink vertices for region overlap
- Returns empty results with warnings for invalid references

## Adding New Tests

To add a new test case to `region_filter_tests.sh`:

1. Define test function following naming convention `test_<feature>`
2. Use helper functions: `run_test`, `print_success`, `print_failure`
3. Add timeout wrappers for povu commands
4. Call test function from `main()`
5. Document expected behavior in comments
