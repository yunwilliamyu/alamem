#!/usr/bin/env bash
THREADS=64
FILELIST="/home/yuny/alamem_benchmark/gtdb_list.txt"
QUERY="1mb.fna"
DATAPOINTS=(1000 2500 5000 10000 20000 30000 40000 50000 60000 70000 80000 85205)
# The last DATAPOINT should be the total length of FILELIST ideally

cd "$(git rev-parse --show-toplevel)" || exit 1
cd scripts

if ! command -v alamem &> /dev/null; then
  echo "Error: alamem is not installed or not on your PATH." >&2
  exit 1
fi

check_gnu_time() {
    # Check 'env time', 'gtime' (common on macOS), or '/usr/bin/time'
    if env time --version 2>&1 | grep -qi "GNU time"; then
        TIME_CMD="env time"
    elif command -v gtime &>/dev/null && gtime --version 2>&1 | grep -qi "GNU time"; then
        TIME_CMD="gtime"
    elif [ -x /usr/bin/time ] && /usr/bin/time --version 2>&1 | grep -qi "GNU time"; then
        TIME_CMD="/usr/bin/time"
    else
        echo "Error: GNU time is required but not installed." >&2
        echo "Install it using 'sudo apt install time' (Linux) or 'brew install time' (macOS)." >&2
        exit 1
    fi
}

check_gnu_time
echo "GNU time found! Running via: $TIME_CMD"

if [ ! -f "$FILELIST" ]; then
  echo "Error: File '$FILELIST' does not exist." >&2
  echo "You should set the parameters at the top of the bash script to point to"
  echo "newline-delimited list of all the FASTAs in your database"
  echo "(e.g. GTDB r214 representative genomes)"
  echo 
  echo "If you need to create this file, you can use"
  echo '  find /dir/to/gtdb/root -type f'
  exit 1
fi

LAST_VAL="${DATAPOINTS[-1]}"
if [ `wc -l < "$FILELIST"` -eq "$LAST_VAL" ]; then
  echo "Match! Last desired benchmark datapoint matches length of filelist."
else
  echo "Mismatch: please edit parameters at top of bash script to ensure that"
  echo "  $FILELIST"
  echo "has the same number of lines as the final benchmark datapoint"
  exit 1
fi


echo Running on subsets of $FILELIST
echo Benchmarking with $THREADS threads
echo Testing runtime and resident memory.

mkdir -p figures



echo -e "N\tTIME\tMEM" > figures/database-size-timing.tsv
for N in "${DATAPOINTS[@]}"
do
  DATABASE="figures/database-$N.txt"
  head -n $N "$FILELIST" > "$DATABASE"
  ALAMEM_CMD="alamem -t $THREADS $DATABASE $QUERY /dev/null" 
  read -r wall_time peak_mem < <($TIME_CMD -f "%e %M" sh -c "$ALAMEM_CMD > /dev/null 2>&1" 2>&1)
  echo -e "$N\t${wall_time}\t${peak_mem}" >> figures/database-size-timing.tsv
  echo -e "$N\t${wall_time}\t${peak_mem}"
done

