download() {
  sra=$1

  if compgen -G "${sra}*" > /dev/null; then
    echo SRA $sra already downloaded...
    return
  fi

  fasterq-dump $sra -q
  pigz ${sra}*.f*q
}

export -f download

parallel -j 4 --line-buffer download ::: $(cat sra.txt)
