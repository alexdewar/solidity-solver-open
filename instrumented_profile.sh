#!/bin/bash

# Instrumented profiling script for Yfd3TET4JOINT with manual timing
# Usage: ./instrumented_profile.sh [num_frames]

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Instrumented Yfd3TET4JOINT Profiler ===${NC}"

# Configuration
NUM_FRAMES=${1:-5}  # Default to 5 frames for profiling
BINARY="./src/bin/Release/DEMPlus"
PROFILE_DIR="instrumented_profile_$(date +%Y%m%d_%H%M%S)"

# Check if binary exists
if [ ! -f "$BINARY" ]; then
    echo -e "${RED}Error: Binary not found at $BINARY${NC}"
    echo "Run: cd src && make release"
    exit 1
fi

# Create profile directory
mkdir -p "$PROFILE_DIR"
echo -e "${GREEN}Profile results will be saved to: $PROFILE_DIR${NC}"

# Copy and modify input files
echo -e "${BLUE}Setting up modified input files for $NUM_FRAMES frames...${NC}"
cp -r examples/* "$PROFILE_DIR/"
sed -i "s/750/$NUM_FRAMES/" "$PROFILE_DIR/layout.txt"
echo -e "${GREEN}Modified simulation to run $NUM_FRAMES frames${NC}"

cd "$PROFILE_DIR"

# Method 1: Instrumented timing run
echo -e "\n${GREEN}=== 1. Instrumented Timing Analysis ===${NC}"
echo "Running instrumented binary with detailed function timing..."

{
    echo "=== Instrumented Run ($NUM_FRAMES frames) ==="
    echo "Start time: $(date)"
    echo "Binary: $BINARY"
    echo "Frames: $NUM_FRAMES"
    echo ""
    
    time ../$BINARY BST.Y3D layout.txt pellet.txt mat.txt 2>&1
    
    echo ""
    echo "End time: $(date)"
} > instrumented_timing.txt 2>&1

# Extract timing reports
echo -e "${BLUE}Extracting detailed timing data...${NC}"

# Extract Yfd3TET4JOINT timing reports
grep -A 10 "Yfd3TET4JOINT Timing Report" instrumented_timing.txt > function_timing_details.txt || echo "No timing reports found"

# Extract overall Yfd performance
grep -E "Yfd:|wall time|frame [0-9]+/" instrumented_timing.txt > overall_yfd_timing.txt

# Method 2: Performance counters with instrumented binary
echo -e "\n${GREEN}=== 2. Performance Counters with Instrumentation ===${NC}"
echo "Collecting hardware counters..."

perf stat -e cycles,instructions,cache-references,cache-misses,branches,branch-misses,page-faults \
  --field-separator=',' --output perf_counters_instrumented.csv \
  ../$BINARY BST.Y3D layout.txt pellet.txt mat.txt > /dev/null 2>&1 || {
  echo -e "${YELLOW}perf stat failed, continuing without detailed counters...${NC}"
}

# Method 3: Thread scaling with instrumented binary
echo -e "\n${GREEN}=== 3. Thread Scaling with Instrumentation ===${NC}"
echo "Testing thread performance with detailed timing..."

cat > thread_scaling_instrumented.txt << EOF
Thread Scaling Analysis with Function-Level Timing
=================================================
EOF

for threads in 1 2 4 8; do
    #if [ $threads -le $(nproc) ]; then
        echo -e "${BLUE}Testing with $threads thread(s)...${NC}"
        export OMP_NUM_THREADS=$threads
        
        echo "=== $threads threads ===" >> thread_scaling_instrumented.txt
        echo "Start: $(date)" >> thread_scaling_instrumented.txt
        
        # Run and capture both overall and function-level timing
        { 
            timeout 300s ../$BINARY BST.Y3D layout.txt pellet.txt mat.txt 2>&1 | \
            grep -E "Yfd:|wall time|Yfd3TET4JOINT Timing Report" | head -20
        } >> thread_scaling_instrumented.txt 2>&1 || {
            echo "Thread test with $threads threads timed out or failed" >> thread_scaling_instrumented.txt
        }
        
        echo "End: $(date)" >> thread_scaling_instrumented.txt
        echo "" >> thread_scaling_instrumented.txt
    #fi
done
unset OMP_NUM_THREADS

# Method 4: Analysis and report generation
echo -e "\n${GREEN}=== 4. Generating Detailed Analysis Report ===${NC}"

cat > detailed_analysis.md << EOF
# Detailed Yfd3TET4JOINT Performance Analysis with Function-Level Timing

## Test Configuration
- Frames analyzed: $NUM_FRAMES
- CPU cores: $(nproc)
- Analysis date: $(date)
- Binary: Instrumented version with manual timing
- Timing precision: Microseconds (gettimeofday)

## Function-Level Timing Results

### Yfd3TET4JOINT Internal Breakdown:
EOF

# Add function timing details if available
if [ -f function_timing_details.txt ] && [ -s function_timing_details.txt ]; then
    echo "\`\`\`" >> detailed_analysis.md
    cat function_timing_details.txt >> detailed_analysis.md
    echo "\`\`\`" >> detailed_analysis.md
else
    echo "No detailed function timing captured - check if timing reports are generated" >> detailed_analysis.md
fi

cat >> detailed_analysis.md << EOF

### Overall Yfd Module Performance:
\`\`\`
$(head -20 overall_yfd_timing.txt)
\`\`\`

## Performance Counter Analysis:
EOF

if [ -f perf_counters_instrumented.csv ]; then
    echo "\`\`\`" >> detailed_analysis.md
    cat perf_counters_instrumented.csv >> detailed_analysis.md
    echo "\`\`\`" >> detailed_analysis.md
    
    # Calculate key metrics
    echo "" >> detailed_analysis.md
    echo "### Key Performance Metrics:" >> detailed_analysis.md
    
    IPC=$(awk -F',' '/instructions/{inst=$1} /cycles/{cyc=$1} END{if(cyc>0) printf "%.3f", inst/cyc}' perf_counters_instrumented.csv)
    CACHE_MISS_RATE=$(awk -F',' '/cache-references/{ref=$1} /cache-misses/{miss=$1} END{if(ref>0) printf "%.2f", (miss/ref)*100}' perf_counters_instrumented.csv)
    
    echo "- Instructions Per Cycle (IPC): $IPC" >> detailed_analysis.md
    echo "- Cache Miss Rate: ${CACHE_MISS_RATE}%" >> detailed_analysis.md
fi

cat >> detailed_analysis.md << EOF

## Thread Scaling Results:
\`\`\`
$(cat thread_scaling_instrumented.txt)
\`\`\`

## Detailed Analysis

### Function Performance Breakdown
Based on the instrumented timing data, we can now identify:

1. **S_N_direction Performance**
   - Time spent in geometry calculations
   - Percentage of total Yfd3TET4JOINT time
   - Per-element average timing

2. **Sigma_tau Performance**  
   - Time spent in stress analysis and integration
   - Percentage of total Yfd3TET4JOINT time
   - Numerical integration overhead

3. **Other Operations**
   - Force accumulation, memory operations, etc.
   - OpenMP overhead and synchronization

### Performance Bottleneck Identification

Based on function-level timing:
EOF

# Add specific analysis based on timing data
if [ -f function_timing_details.txt ] && [ -s function_timing_details.txt ]; then
    # Try to extract percentages if available
    S_N_PERCENT=$(grep "S_N_direction" function_timing_details.txt | tail -1 | grep -o '([0-9.]*%)' | tr -d '()')
    SIGMA_PERCENT=$(grep "Sigma_tau" function_timing_details.txt | tail -1 | grep -o '([0-9.]*%)' | tr -d '()')
    
    cat >> detailed_analysis.md << EOF

- **Primary bottleneck**: $(if [[ "${SIGMA_PERCENT%\%}" > "${S_N_PERCENT%\%}" ]] 2>/dev/null; then echo "Sigma_tau ($SIGMA_PERCENT)"; else echo "S_N_direction ($S_N_PERCENT)"; fi)
- **Secondary bottleneck**: $(if [[ "${SIGMA_PERCENT%\%}" > "${S_N_PERCENT%\%}" ]] 2>/dev/null; then echo "S_N_direction ($S_N_PERCENT)"; else echo "Sigma_tau ($SIGMA_PERCENT)"; fi)
EOF
fi

cat >> detailed_analysis.md << EOF

### Optimization Recommendations

1. **Immediate Actions (if Sigma_tau dominates)**:
   - Vectorize numerical integration loops
   - Precompute expensive transcendental functions
   - Optimize Mohr-Coulomb failure criterion evaluation

2. **Immediate Actions (if S_N_direction dominates)**:
   - Cache geometric calculations
   - Optimize coordinate transformations
   - Reduce redundant sqrt/division operations

3. **Memory Access Optimization**:
   - Improve cache locality for element traversal
   - Consider data structure reorganization
   - Prefetch memory for predictable access patterns

4. **Thread Scaling Improvements**:
   - Reduce critical section overhead
   - Better load balancing for non-uniform element processing
   - Consider lock-free data structures

### Next Steps for Further Optimization

1. **Profile individual mathematical operations**
2. **Analyze memory access patterns with cachegrind**
3. **Test SIMD implementations for hot loops**
4. **Evaluate GPU acceleration potential**

## Quick Commands for Analysis:

\`\`\`bash
# View function-level timing breakdown
cat function_timing_details.txt

# Compare thread scaling
grep -A 2 "=== [0-9]* threads ===" thread_scaling_instrumented.txt

# Calculate performance efficiency
grep "real" thread_scaling_instrumented.txt

# View cache performance
cat perf_counters_instrumented.csv
\`\`\`
EOF

# Generate summary
echo -e "\n${GREEN}=== Performance Summary ===${NC}"

echo -e "${YELLOW}Function-Level Timing (Latest Report):${NC}"
if [ -f function_timing_details.txt ] && [ -s function_timing_details.txt ]; then
    tail -10 function_timing_details.txt
else
    echo "No detailed timing data captured"
fi

echo -e "\n${YELLOW}Overall Performance:${NC}"
if [ -f overall_yfd_timing.txt ]; then
    grep "Yfd:" overall_yfd_timing.txt | head -5
fi

echo -e "\n${YELLOW}Thread Scaling Summary:${NC}"
if [ -f thread_scaling_instrumented.txt ]; then
    grep -E "=== [0-9]+ threads ===|real" thread_scaling_instrumented.txt | head -10
fi

echo -e "\n${GREEN}=== Analysis Complete ===${NC}"
echo -e "Results saved in: ${YELLOW}$(pwd)${NC}"
echo -e "\nKey files to review:"
echo -e "  ${GREEN}detailed_analysis.md${NC}          - Comprehensive analysis"
echo -e "  ${GREEN}function_timing_details.txt${NC}   - Function-level timing"
echo -e "  ${GREEN}instrumented_timing.txt${NC}       - Full execution log"
echo -e "  ${GREEN}thread_scaling_instrumented.txt${NC} - Threading performance"

echo -e "\n${GREEN}Instrumented profiling complete!${NC}"
echo -e "${BLUE}Review detailed_analysis.md for specific optimization recommendations${NC}"
