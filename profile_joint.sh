#!/bin/bash

# Comprehensive perf profiling script for Yfd3TET4JOINT function
# Usage: ./profile_joint.sh [input_files...]

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Solidity Yfd3TET4JOINT Profiler ===${NC}"

# Check if perf is installed
if ! command -v perf &> /dev/null; then
    echo -e "${RED}Error: perf is not installed. Install with: sudo apt-get install linux-tools-common linux-tools-generic${NC}"
    exit 1
fi

# Check if binary exists
BINARY="./src/bin/Release/DEMPlus"
if [ ! -f "$BINARY" ]; then
    echo -e "${RED}Error: Binary not found at $BINARY${NC}"
    echo "Run: cd src && make release"
    exit 1
fi

# Create profile output directory
PROFILE_DIR="profile_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$PROFILE_DIR"
echo -e "${GREEN}Profile results will be saved to: $PROFILE_DIR${NC}"

# Set up input files (use examples if none provided)
if [ $# -eq 0 ]; then
    echo -e "${YELLOW}No input files specified, using examples...${NC}"
    cd examples
    INPUT_FILES="BST.Y3D layout.txt pellet.txt mat.txt"
    WORK_DIR="$(pwd)"
else
    INPUT_FILES="$@"
    WORK_DIR="$(pwd)"
fi

echo -e "${GREEN}Input files: $INPUT_FILES${NC}"

# Function to run perf with specific options
run_perf_analysis() {
    local name="$1"
    local perf_opts="$2"
    local output_file="$PROFILE_DIR/${name}.txt"
    
    echo -e "\n${YELLOW}Running $name analysis...${NC}"
    
    # Run perf with timeout to prevent hanging
    timeout 300s perf $perf_opts -- $BINARY $INPUT_FILES > "$output_file" 2>&1 || {
        echo -e "${RED}Warning: $name analysis timed out or failed${NC}"
        return 1
    }
    
    echo -e "${GREEN}✓ $name results saved to $output_file${NC}"
}

# 1. CPU Usage Profiling - Focus on Yfd3TET4JOINT
echo -e "\n${GREEN}=== 1. CPU Usage Analysis ===${NC}"
run_perf_analysis "cpu_profile" \
    "record -g -F 999 --call-graph=dwarf -e cycles -o $PROFILE_DIR/cpu.data"

if [ -f "$PROFILE_DIR/cpu.data" ]; then
    perf report -i "$PROFILE_DIR/cpu.data" --stdio > "$PROFILE_DIR/cpu_report.txt"
    perf report -i "$PROFILE_DIR/cpu.data" --stdio --sort=symbol | grep -A 20 -B 5 "Yfd3TET4JOINT" > "$PROFILE_DIR/joint_function_focus.txt" || echo "Function not found in top results"
fi

# 2. Cache Analysis
echo -e "\n${GREEN}=== 2. Cache Performance Analysis ===${NC}"
run_perf_analysis "cache_profile" \
    "stat -e cache-references,cache-misses,LLC-loads,LLC-load-misses,L1-dcache-loads,L1-dcache-load-misses,L1-icache-load-misses -o $PROFILE_DIR/cache_stats.txt"

# 3. Memory Bandwidth Analysis
echo -e "\n${GREEN}=== 3. Memory Bandwidth Analysis ===${NC}"
run_perf_analysis "memory_profile" \
    "stat -e cycles,instructions,cache-references,cache-misses,bus-cycles,mem_load_retired.l3_miss,mem_load_retired.l3_hit -o $PROFILE_DIR/memory_stats.txt"

# 4. Branch Prediction Analysis
echo -e "\n${GREEN}=== 4. Branch Prediction Analysis ===${NC}"
run_perf_analysis "branch_profile" \
    "stat -e branches,branch-misses,branch-load-misses -o $PROFILE_DIR/branch_stats.txt"

# 5. OpenMP/Threading Analysis
echo -e "\n${GREEN}=== 5. Threading Analysis ===${NC}"
run_perf_analysis "thread_profile" \
    "record -e sched:sched_switch,sched:sched_wakeup -a -g -o $PROFILE_DIR/thread.data"

if [ -f "$PROFILE_DIR/thread.data" ]; then
    perf script -i "$PROFILE_DIR/thread.data" > "$PROFILE_DIR/thread_trace.txt"
fi

# 6. Detailed Function-level Profiling with Annotations
echo -e "\n${GREEN}=== 6. Function Annotation Analysis ===${NC}"
if command -v perf &> /dev/null && [ -f "$PROFILE_DIR/cpu.data" ]; then
    echo "Generating annotated source for Yfd3TET4JOINT..."
    perf annotate -i "$PROFILE_DIR/cpu.data" --stdio Yfd3TET4JOINT > "$PROFILE_DIR/yfd3tet4joint_annotated.txt" 2>/dev/null || {
        echo -e "${YELLOW}Could not generate function annotation (function may not be in top hotspots)${NC}"
    }
fi

# 7. Top Functions Report
echo -e "\n${GREEN}=== 7. Generating Summary Reports ===${NC}"

cat > "$PROFILE_DIR/README.md" << EOF
# Yfd3TET4JOINT Performance Profile Results

Generated: $(date)

## Files Description:

### CPU Profiling:
- **cpu_report.txt**: Complete CPU usage breakdown
- **joint_function_focus.txt**: Specific focus on Yfd3TET4JOINT function
- **yfd3tet4joint_annotated.txt**: Source-level annotations (if available)

### Memory Profiling:
- **cache_stats.txt**: L1/L2/L3 cache performance
- **memory_stats.txt**: Memory bandwidth and access patterns

### Other Analysis:
- **branch_stats.txt**: Branch prediction performance
- **thread_trace.txt**: Thread scheduling and synchronization
- **summary_report.txt**: Key findings and recommendations

## Key Metrics to Check:

1. **CPU Hotspots**: Look for Yfd3TET4JOINT in cpu_report.txt
2. **Cache Misses**: High cache miss rates indicate memory bottlenecks
3. **Branch Mispredictions**: High rates suggest unpredictable code paths
4. **Thread Contention**: Check thread_trace.txt for synchronization issues

## Quick Commands:

\`\`\`bash
# View top functions consuming CPU
head -50 cpu_report.txt

# Check cache miss rate
grep -E "(cache-misses|cache-references)" cache_stats.txt

# Find OpenMP critical section overhead
grep -i "critical\\|atomic\\|omp" cpu_report.txt
\`\`\`
EOF

# Generate a summary report
cat > "$PROFILE_DIR/summary_report.txt" << EOF
# Performance Analysis Summary for Yfd3TET4JOINT

## Analysis Overview:
- Binary: $BINARY
- Input: $INPUT_FILES
- Analysis Date: $(date)

## Key Areas to Investigate:

1. **Function Hotspots**
   - Check if Yfd3TET4JOINT appears in top CPU consumers
   - Look for S_N_direction and Sigma_tau subfunctions

2. **Memory Access Patterns**
   - Array accesses to d2ejfjx, d2ejfjy, d2ejfjz
   - Node coordinate lookups (d1nccx, d1nccy, d1nccz)
   - Element connectivity traversal (i2elto)

3. **OpenMP Overhead**
   - Critical section contention when elements fail
   - Load balancing across threads
   - False sharing in shared arrays

4. **Computational Bottlenecks**
   - Mathematical operations (sqrt, exp, pow)
   - Integration loop (3 points per element)
   - Failure criterion evaluation

## Next Steps:
1. Review cpu_report.txt for function rankings
2. Check cache_stats.txt for memory efficiency  
3. Examine branch_stats.txt for prediction issues
4. Consider algorithmic improvements based on findings
EOF

echo -e "\n${GREEN}=== Analysis Complete ===${NC}"
echo -e "Results saved in: ${YELLOW}$PROFILE_DIR${NC}"
echo -e "\nTo view results:"
echo -e "  ${GREEN}cd $PROFILE_DIR${NC}"
echo -e "  ${GREEN}cat README.md${NC}"
echo -e "  ${GREEN}less cpu_report.txt${NC}"

# Quick summary of file sizes (indicates if data was captured)
echo -e "\n${YELLOW}Generated Files:${NC}"
ls -la "$PROFILE_DIR"/ | grep -E '\.(txt|data)$' | awk '{printf "  %-30s %8s\n", $9, $5}'

echo -e "\n${GREEN}Profile analysis complete!${NC}"