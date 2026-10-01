#!/usr/bin/env bash
cd "$(git rev-parse --show-toplevel)" || exit 1
cd scripts
echo "Compiling ani_simulator in release mode"
RUSTFLAGS="-C target-cpu=native" cargo build --example ani_simulator --release

mkdir -p figures

echo 'Testing accuracy for ANI in {90..100} and fragment length in {50..2000..50} and k in {11,13,15}'

for k in 11 13 15
do

  echo -e "ANI\tLength\tMapped Proportion\tMean Absolute Error\tMean Bias\tStandard Error\tMean Length" > figures/accuracy-$k.tsv
  for ani in {90..100}
  do
    for l in {50..2000..50}
    do
      (
      LINE=`../target/release/examples/ani_simulator ./1mb.fna 10000 $l $ani -k $k -m 1000 --min-ani 90 | tail -n 5 | cut -f 2 -d ':' | tr -d '[:blank:]' | paste -sd "\t" | sed 's/reads//'`
      echo -e "$ani\t$l\t$LINE" >> figures/accuracy-$k.tsv
      echo "Simulated ANI=$ani and Length=$l and k=$k"
      ) &
    done
    wait
  done
  { head -n 1 figures/accuracy-$k.tsv; tail -n +2 figures/accuracy-$k.tsv | sort -k1,1n -k2,2n; } > figures/tmp.tsv && mv figures/tmp.tsv figures/accuracy-$k.tsv

  if command -v python > /dev/null 2>&1; then
    cd figures
    echo "Generating figure using Python"
    python ../accuracy-bench-plot.py accuracy-$k.tsv
    cd ..
  else
    echo "Python 3 is not installed. Cannot generate figure."
  fi

done

