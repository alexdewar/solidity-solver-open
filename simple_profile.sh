#!/bin/bash

# Simple profiling script for Yfd3TET4JOINT that works with restricted perf access
# Usage: ./simple_profile.sh [input_files...]

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Simple Yfd3TET4JOINT Profiler ===${NC}"

# Check if binary exists
BINARY="./src/bin/Release/DEMPlus"
if [ ! -f "$BINARY" ]; then
    echo -e "${RED}Error: Binary not found at $BINARY${NC}"
    echo "Run: cd src && make release"
    exit 1
fi

# Create profile output directory
PROFILE_DIR="simple_profile_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$PROFILE_DIR"
echo -e "${GREEN}Profile results will be saved to: $PROFILE_DIR${NC}"

# Set up input files (use examples if none provided)
if [ $# -eq 0 ]; then
    echo -e "${YELLOW}No input files specified, using examples...${NC}"
    cd examples
    INPUT_FILES="BST.Y3D layout.txt pellet.txt mat.txt"
    WORK_DIR="$(pwd)"
    BINARY="../$BINARY"
else
    INPUT_FILES="$@"
    WORK_DIR="$(pwd)"
fi

echo -e "${GREEN}Input files: $INPUT_FILES${NC}"

# Method 1: Time-based analysis with gprof (if available)
echo -e "\n${GREEN}=== 1. Timing Analysis ===${NC}"

# First, let's just time the execution
echo "Running basic timing test..."
time $BINARY $INPUT_FILES > "$PROFILE_DIR/execution_output.txt" 2> "$PROFILE_DIR/timing_basic.txt"

# Method 2: Simple perf stat (works with paranoid=2)
echo -e "\n${GREEN}=== 2. Basic Performance Counters ===${NC}"
echo "Collecting basic performance statistics..."

perf stat -e cycles,instructions,cache-references,cache-misses,branches,branch-misses \
  $BINARY $INPUT_FILES > /dev/null 2> "$PROFILE_DIR/perf_stat.txt" || {
  echo -e "${YELLOW}perf stat failed, continuing without it...${NC}"
}

# Method 3: Valgrind callgrind (if available)
if command -v valgrind &> /dev/null; then
    echo -e "\n${GREEN}=== 3. Valgrind Callgrind Analysis ===${NC}"
    echo "Running callgrind analysis (this may take a while)..."
    
    timeout 300s valgrind --tool=callgrind --callgrind-out-file="$PROFILE_DIR/callgrind.out" \
      $BINARY $INPUT_FILES > /dev/null 2> "$PROFILE_DIR/valgrind_log.txt" || {
      echo -e "${YELLOW}Valgrind analysis timed out or failed${NC}"
    }
    
    if [ -f "$PROFILE_DIR/callgrind.out" ]; then
        echo "Generating callgrind report..."
        callgrind_annotate "$PROFILE_DIR/callgrind.out" > "$PROFILE_DIR/callgrind_report.txt" 2>/dev/null || {
            echo -e "${YELLOW}callgrind_annotate not available${NC}"
        }
    fi
else
    echo -e "${YELLOW}Valgrind not available, skipping callgrind analysis${NC}"
fi

# Method 4: Manual code instrumentation suggestion
echo -e "\n${GREEN}=== 4. Creating Manual Instrumentation Guide ===${NC}"

cat > "$PROFILE_DIR/manual_instrumentation.c" << 'EOF'
/* Manual timing instrumentation for Yfd3TET4JOINT */
/* Add these includes at the top of Y3Dfd.c */

#include <sys/time.h>
#include <stdio.h>

/* Add these variables at the top of Yfd3TET4JOINT function */
struct timeval start_time, end_time;
double total_time = 0.0;
static int call_count = 0;
double s_n_direction_time = 0.0;
double sigma_tau_time = 0.0;
double force_accumulation_time = 0.0;

/* At the beginning of Yfd3TET4JOINT function */
gettimeofday(&start_time, NULL);
call_count++;

/* Before S_N_direction call */
struct timeval s_n_start;
gettimeofday(&s_n_start, NULL);

/* After S_N_direction call */
struct timeval s_n_end;
gettimeofday(&s_n_end, NULL);
s_n_direction_time += (s_n_end.tv_sec - s_n_start.tv_sec) + 
                      (s_n_end.tv_usec - s_n_start.tv_usec) / 1000000.0;

/* Before Sigma_tau call */
struct timeval sigma_start;
gettimeofday(&sigma_start, NULL);

/* After Sigma_tau call */
struct timeval sigma_end;
gettimeofday(&sigma_end, NULL);
sigma_tau_time += (sigma_end.tv_sec - sigma_start.tv_sec) + 
                  (sigma_end.tv_usec - sigma_start.tv_usec) / 1000000.0;

/* Before force accumulation loop (line ~1320) */
struct timeval force_start;
gettimeofday(&force_start, NULL);

/* After force accumulation loop */
struct timeval force_end;
gettimeofday(&force_end, NULL);
force_accumulation_time += (force_end.tv_sec - force_start.tv_sec) + 
                           (force_end.tv_usec - force_start.tv_usec) / 1000000.0;

/* At the end of Yfd3TET4JOINT function */
gettimeofday(&end_time, NULL);
total_time += (end_time.tv_sec - start_time.tv_sec) + 
              (end_time.tv_usec - start_time.tv_usec) / 1000000.0;

/* Print timing info every 100 calls */
if (call_count % 100 == 0) {
    printf("Yfd3TET4JOINT call %d: total=%.6fs, S_N_direction=%.6fs, Sigma_tau=%.6fs, force_accum=%.6fs\n",
           call_count, total_time, s_n_direction_time, sigma_tau_time, force_accumulation_time);
}
EOF

# Method 5: Simple threading analysis
echo -e "\n${GREEN}=== 5. Threading Analysis ===${NC}"

# Check OpenMP environment
echo "OpenMP Settings:" > "$PROFILE_DIR/openmp_analysis.txt"
echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS:-not set}" >> "$PROFILE_DIR/openmp_analysis.txt"
echo "Available CPU cores: $(nproc)" >> "$PROFILE_DIR/openmp_analysis.txt"

# Run with different thread counts if OMP is available
if [ -z "$OMP_NUM_THREADS" ]; then
    echo "Testing different thread counts..."
    
    for threads in 1 2 4; do
        if [ $threads -le $(nproc) ]; then
            echo "Testing with $threads thread(s)..."
            export OMP_NUM_THREADS=$threads
            echo "=== $threads threads ===" >> "$PROFILE_DIR/openmp_analysis.txt"
            time $BINARY $INPUT_FILES > /dev/null 2>> "$PROFILE_DIR/openmp_analysis.txt"
        fi
    done
    unset OMP_NUM_THREADS
fi

# Create analysis report
cat > "$PROFILE_DIR/analysis_report.md" << EOF
# Simple Performance Analysis Report

## Analysis Date
$(date)

## Binary Information
- Executable: $BINARY
- Input Files: $INPUT_FILES
- Compiler Flags: -O3 -g -fopenmp -funroll-loops -ftree-vectorize

## Available Results

### 1. Basic Timing (timing_basic.txt)
Contains wall-clock time, user time, and system time for the complete execution.

### 2. Performance Counters (perf_stat.txt)
Basic hardware counters including:
- CPU cycles and instructions
- Cache references and misses
- Branch predictions and misses

### 3. Callgrind Analysis (if available)
- callgrind.out: Raw callgrind data
- callgrind_report.txt: Function call costs

### 4. Threading Analysis (openmp_analysis.txt)
Comparison of execution times with different thread counts.

## Quick Analysis Commands

\`\`\`bash
# Check basic execution time
cat timing_basic.txt

# View performance counters
cat perf_stat.txt

# Check cache miss rate
grep -E "cache.*(references|misses)" perf_stat.txt

# View function costs (if callgrind available)
head -50 callgrind_report.txt | grep -E "(Yfd3TET4JOINT|S_N_direction|Sigma_tau)"

# Compare thread performance
grep "real" openmp_analysis.txt
\`\`\`

## Key Metrics to Look For

1. **High Cache Miss Rate**: > 10% indicates memory bottlenecks
2. **Branch Misprediction**: > 5% suggests unpredictable code paths  
3. **Low IPC (Instructions Per Cycle)**: < 1.0 indicates compute bottlenecks
4. **Thread Scaling**: Speedup should be close to number of threads

## Recommendations for Further Analysis

1. **Add Manual Timing**: Use the code in manual_instrumentation.c
2. **Profile with Root Access**: Use perf record for detailed analysis
3. **Memory Access Patterns**: Use valgrind --tool=massif
4. **Code Coverage**: Identify hot loops within Yfd3TET4JOINT

## Expected Bottlenecks in Yfd3TET4JOINT

Based on code analysis, likely bottlenecks:
1. S_N_direction geometry calculations (sqrt, divisions)
2. Sigma_tau numerical integration (3 points per element)
3. Critical sections when elements fail
4. Force accumulation loop (sequential)
5. Memory access patterns for large meshes
EOF

echo -e "\n${GREEN}=== Analysis Complete ===${NC}"
echo -e "Results saved in: ${YELLOW}$PROFILE_DIR${NC}"
echo -e "\nTo view results:"
echo -e "  ${GREEN}cd $PROFILE_DIR${NC}"
echo -e "  ${GREEN}cat analysis_report.md${NC}"
echo -e "  ${GREEN}cat timing_basic.txt${NC}"

# Show quick summary
echo -e "\n${YELLOW}Quick Results Summary:${NC}"
if [ -f "$PROFILE_DIR/timing_basic.txt" ]; then
    echo "Execution Time:"
    cat "$PROFILE_DIR/timing_basic.txt" | tail -3
fi

if [ -f "$PROFILE_DIR/perf_stat.txt" ]; then
    echo -e "\nKey Performance Counters:"
    grep -E "(cycles|instructions|cache-misses|branch-misses)" "$PROFILE_DIR/perf_stat.txt" | head -6
fi

echo -e "\n${GREEN}Simple profiling complete! Check $PROFILE_DIR for detailed results.${NC}"