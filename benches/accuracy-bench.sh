#!/usr/bin/env bash
cd "$(git rev-parse --show-toplevel)" || exit 1
cd benches
echo "Compiling ani_simulator in release mode"
RUSTFLAGS="-C target-cpu=native" cargo build --example ani_simulator --release

mkdir -p figures

echo 'Testing accuracy for ANI in {90..100} and fragment length in {50..2000..50}'
echo -e "ANI\tLength\tMapped Proportion\tMean Absolute Error\tMean Bias\tStandard Error\tMean Length" > figures/accuracy.tsv
for ani in {90..100}
do
  for l in {50..2000..50}
  do
    LINE=`../target/release/examples/ani_simulator ./1mb.fna 1000 $l $ani -k 11 -m 1000 | tail -n 5 | cut -f 2 -d ':' | tr -d '[:blank:]' | paste -sd "\t" | sed 's/reads//'`
    echo -e "$ani\t$l\t$LINE" >> figures/accuracy.tsv
    echo "Simulated ANI=$ani and Length=$l"
  done
done

if command -v python3 > /dev/null 2>&1; then
  cd figures
  echo "Generating figure using Python"
  python ../accuracy-bench-plot.py
else
  echo "Python 3 is not installed. Cannot generate figure."
fi


