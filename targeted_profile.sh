#!/bin/bash

# Targeted profiling script for Yfd3TET4JOINT based on working binary
# Usage: ./targeted_profile.sh [num_frames]

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Targeted Yfd3TET4JOINT Profiler ===${NC}"

# Configuration
NUM_FRAMES=${1:-10}  # Default to 10 frames for faster profiling
BINARY="./src/bin/Release/DEMPlus"
PROFILE_DIR="targeted_profile_$(date +%Y%m%d_%H%M%S)"

# Check if binary exists
if [ ! -f "$BINARY" ]; then
    echo -e "${RED}Error: Binary not found at $BINARY${NC}"
    echo "Run: cd src && make release"
    exit 1
fi

# Create profile directory
mkdir -p "$PROFILE_DIR"
echo -e "${GREEN}Profile results will be saved to: $PROFILE_DIR${NC}"

# Copy input files for modification
echo -e "${BLUE}Setting up modified input files for $NUM_FRAMES frames...${NC}"
cp -r examples/* "$PROFILE_DIR/"

# Modify layout.txt to reduce simulation time for profiling
sed -i "s/750/$NUM_FRAMES/" "$PROFILE_DIR/layout.txt"
echo -e "${GREEN}Modified simulation to run $NUM_FRAMES frames instead of 750${NC}"

cd "$PROFILE_DIR"

# Method 1: Basic timing with detailed output
echo -e "\n${GREEN}=== 1. Basic Performance Measurement ===${NC}"
echo "Running baseline timing with $NUM_FRAMES frames..."

# Capture detailed timing output
{
    echo "=== Baseline Run ($NUM_FRAMES frames) ==="
    echo "Start time: $(date)"
    time ../$BINARY BST.Y3D layout.txt pellet.txt mat.txt 2>&1
    echo "End time: $(date)"
} > baseline_timing.txt 2>&1

# Extract Yfd timing information
echo -e "${BLUE}Extracting Yfd performance data...${NC}"
grep -E "Yfd:|wall time|frame [0-9]+/" baseline_timing.txt > yfd_timing_analysis.txt

# Method 2: perf stat with basic counters (works without root)
echo -e "\n${GREEN}=== 2. Hardware Performance Counters ===${NC}"
echo "Collecting performance statistics..."

perf stat -e cycles,instructions,cache-references,cache-misses,branches,branch-misses,page-faults \
  --field-separator=',' --output perf_counters.csv \
  ../$BINARY BST.Y3D layout.txt pellet.txt mat.txt > /dev/null 2>&1 || {
  echo -e "${YELLOW}perf stat failed, trying simplified version...${NC}"
  
  # Fallback to simpler counters
  perf stat -e cycles,instructions \
    --field-separator=',' --output perf_counters_simple.csv \
    ../$BINARY BST.Y3D layout.txt pellet.txt mat.txt > /dev/null 2>&1 || {
    echo -e "${RED}perf unavailable${NC}"
  }
}

# Method 3: Thread scaling analysis
echo -e "\n${GREEN}=== 3. Thread Scaling Analysis ===${NC}"
echo "Testing thread scalability..."

cat > thread_scaling.txt << EOF
Thread Scaling Analysis for Yfd3TET4JOINT
==========================================
EOF

for threads in 1 2 4 8 12; do
    if [ $threads -le $(nproc) ]; then
        echo -e "${BLUE}Testing with $threads thread(s)...${NC}"
        export OMP_NUM_THREADS=$threads
        
        echo "=== $threads threads ===" >> thread_scaling.txt
        echo "Start: $(date)" >> thread_scaling.txt
        
        # Run and capture timing
        { time ../$BINARY BST.Y3D layout.txt pellet.txt mat.txt 2>&1 | \
          grep -E "Yfd:|wall time" | head -10; } >> thread_scaling.txt 2>&1
        
        echo "End: $(date)" >> thread_scaling.txt
        echo "" >> thread_scaling.txt
    fi
done
unset OMP_NUM_THREADS

# Method 4: Memory usage analysis with /proc/meminfo
echo -e "\n${GREEN}=== 4. Memory Usage Analysis ===${NC}"
echo "Monitoring memory usage during execution..."

# Background memory monitor
(
    echo "Memory usage over time:" > memory_usage.txt
    while true; do
        echo "$(date '+%H:%M:%S'): $(grep -E 'MemAvailable|MemFree' /proc/meminfo | tr '\n' ' ')" >> memory_usage.txt
        sleep 1
    done
) &
MONITOR_PID=$!

# Run the simulation
../$BINARY BST.Y3D layout.txt pellet.txt mat.txt > memory_run_output.txt 2>&1

# Stop memory monitoring
kill $MONITOR_PID 2>/dev/null || true

# Method 5: Create analysis summary
echo -e "\n${GREEN}=== 5. Generating Analysis Report ===${NC}"

cat > analysis_summary.md << EOF
# Targeted Performance Analysis for Yfd3TET4JOINT

## Test Configuration
- Frames analyzed: $NUM_FRAMES (reduced from 750 for faster profiling)
- CPU cores available: $(nproc)
- Analysis date: $(date)
- Binary: $BINARY

## Key Findings

### Timing Analysis
Based on baseline run with $NUM_FRAMES frames:

EOF

# Extract key timing metrics
if [ -f yfd_timing_analysis.txt ]; then
    echo "#### Yfd Function Performance:" >> analysis_summary.md
    echo "\`\`\`" >> analysis_summary.md
    head -20 yfd_timing_analysis.txt >> analysis_summary.md
    echo "\`\`\`" >> analysis_summary.md
    echo "" >> analysis_summary.md
fi

# Analyze performance counters
if [ -f perf_counters.csv ]; then
    echo "#### Hardware Performance Counters:" >> analysis_summary.md
    echo "\`\`\`" >> analysis_summary.md
    cat perf_counters.csv >> analysis_summary.md
    echo "\`\`\`" >> analysis_summary.md
    echo "" >> analysis_summary.md
fi

cat >> analysis_summary.md << EOF

### Thread Scaling Results
The following shows how Yfd performance scales with different thread counts:

\`\`\`
$(grep -A 2 "=== [0-9]* threads ===" thread_scaling.txt | head -20)
\`\`\`

## Performance Bottleneck Analysis

### Expected Bottlenecks in Yfd3TET4JOINT:

1. **Geometry Calculations (S_N_direction)**
   - Heavy use of sqrt(), divisions
   - Coordinate transformations
   - Normal vector calculations

2. **Stress Analysis (Sigma_tau)**
   - Numerical integration (3 points per element)
   - Mohr-Coulomb failure criterion evaluation
   - Exponential and power functions

3. **Memory Access Patterns**
   - Large arrays: d2ejfjx, d2ejfjy, d2ejfjz
   - Element connectivity traversal
   - Node coordinate lookups

4. **Synchronization Overhead**
   - Critical sections when elements fail
   - Atomic operations on shared data
   - Thread load imbalancing

### Recommended Next Steps:

1. **Profile Individual Functions**
   - Add timing instrumentation to S_N_direction
   - Measure Sigma_tau computation time
   - Profile force accumulation loop

2. **Memory Optimization**
   - Analyze cache miss patterns
   - Consider data structure reorganization
   - Evaluate memory access locality

3. **Algorithm Optimization**
   - Vectorize mathematical operations
   - Reduce redundant calculations
   - Optimize critical sections

4. **Scaling Improvements**
   - Better load balancing strategies
   - Reduce synchronization overhead
   - Consider GPU acceleration for compute-heavy parts

## Quick Commands for Further Analysis:

\`\`\`bash
# View detailed timing breakdown
grep "Yfd:" baseline_timing.txt

# Check thread scaling efficiency
grep "real" thread_scaling.txt

# Monitor memory usage pattern
tail -20 memory_usage.txt

# Calculate cache miss rate (if available)
awk -F',' '/cache-references/{ref=\$1} /cache-misses/{miss=\$1} END{if(ref>0) print "Cache miss rate:", (miss/ref)*100 "%"}' perf_counters.csv
\`\`\`

EOF

# Generate quick statistics
echo -e "\n${GREEN}=== Quick Results Summary ===${NC}"

echo -e "${YELLOW}Yfd Performance (from baseline run):${NC}"
if [ -f yfd_timing_analysis.txt ]; then
    grep "Yfd:" yfd_timing_analysis.txt | head -5
fi

echo -e "\n${YELLOW}Thread Scaling Summary:${NC}"
if [ -f thread_scaling.txt ]; then
    grep -B 1 -A 1 "real" thread_scaling.txt | head -10
fi

echo -e "\n${YELLOW}Performance Counter Summary:${NC}"
if [ -f perf_counters.csv ]; then
    echo "Instructions per cycle (IPC):"
    awk -F',' '/instructions/{inst=$1} /cycles/{cyc=$1} END{if(cyc>0) printf "%.2f\n", inst/cyc}' perf_counters.csv
    
    echo "Cache miss rate:"
    awk -F',' '/cache-references/{ref=$1} /cache-misses/{miss=$1} END{if(ref>0) printf "%.2f%%\n", (miss/ref)*100}' perf_counters.csv
fi

echo -e "\n${GREEN}=== Analysis Complete ===${NC}"
echo -e "Results saved in: ${YELLOW}$(pwd)${NC}"
echo -e "\nTo view results:"
echo -e "  ${GREEN}cat analysis_summary.md${NC}"
echo -e "  ${GREEN}less baseline_timing.txt${NC}"
echo -e "  ${GREEN}cat yfd_timing_analysis.txt${NC}"

echo -e "\n${GREEN}Targeted profiling complete!${NC}"
echo -e "${BLUE}Next: Review analysis_summary.md for optimization opportunities${NC}"