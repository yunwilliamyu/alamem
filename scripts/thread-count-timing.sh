#!/usr/bin/env bash
THREADS=(1 2 4 8 16 32 64)
FILELIST="/home/yuny/alamem_benchmark/gtdb_list.txt"
QUERY="1mb.fna"
DATAPOINTS=10000 # Running thread benchmarking is slow, so we subsample the database

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


echo Running on subsets of $FILELIST
echo Benchmarking with $THREADS threads
echo Testing runtime and resident memory.

mkdir -p figures



echo -e "THREADS\tTIME\tMEM" > figures/thread-count-timing.tsv
for T in "${THREADS[@]}"
do
  DATABASE="figures/database-$DATAPOINTS.txt"
  head -n $DATAPOINTS "$FILELIST" > "$DATABASE"
  ALAMEM_CMD="alamem -t $T $DATABASE $QUERY /dev/null" 
  read -r wall_time peak_mem < <($TIME_CMD -f "%e %M" sh -c "$ALAMEM_CMD > /dev/null 2>&1" 2>&1)
  echo -e "$T\t${wall_time}\t${peak_mem}" >> figures/thread-count-timing.tsv
  echo -e "$T\t${wall_time}\t${peak_mem}"
  rm $DATABASE
done

