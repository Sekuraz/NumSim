BEGIN {
  FS="\t";
  OFS="\t";
  N=50;
}
function abs(x) { return (x<0) ? -x : x }
{
  re+=$1; reSq+=$1^2;
  u1+=$4; u1Sq+=$4^2; v1+=$5; v1Sq+=$5^2;
  u2+=$6; u2Sq+=$6^2; v2+=$7; v2Sq+=$7^2;
  u3+=$8; u3Sq+=$8^2; v3+=$9; v3Sq+=$9^2;
}
NR%N == 0 {
  i = int( NR/N );
  reM[i]=re/NR;
  u1M[i]=u1/NR;
  v1M[i]=v1/NR;
  u2M[i]=u2/NR;
  v2M[i]=v2/NR;
  u3M[i]=u3/NR;
  v3M[i]=v3/NR;
  reSD[i]=sqrt((reSq - re^2/NR)/(NR-1));
  u1SD[i]=sqrt((u1Sq - u1^2/NR)/(NR-1));
  v1SD[i]=sqrt((v1Sq - v1^2/NR)/(NR-1));
  u2SD[i]=sqrt((u2Sq - u2^2/NR)/(NR-1));
  v2SD[i]=sqrt((v2Sq - v2^2/NR)/(NR-1));
  u3SD[i]=sqrt((u3Sq - u3^2/NR)/(NR-1));
  v3SD[i]=sqrt((v3Sq - v3^2/NR)/(NR-1));
}
END {
  Nend = int( NR/N );
  for(i=1; i<Nend; i++) {
    print i*N,
          abs(reM[i]-1500)/abs(reM[1]-1500), abs(reSD[i]-166.666666666667)/abs(reSD[1]-166.666666666667),
          abs(u1M[i]-u1M[Nend])/abs(u1M[1]-u1M[Nend]), abs(u1SD[i]-u1SD[Nend])/abs(u1SD[1]-u1SD[Nend]),
          abs(v1M[i]-v1M[Nend])/abs(v1M[1]-v1M[Nend]), abs(v1SD[i]-v1SD[Nend])/abs(v1SD[1]-v1SD[Nend]),
          abs(u2M[i]-u2M[Nend])/abs(u2M[1]-u2M[Nend]), abs(u2SD[i]-u2SD[Nend])/abs(u2SD[1]-u2SD[Nend]),
          abs(v2M[i]-v2M[Nend])/abs(v2M[1]-v2M[Nend]), abs(v2SD[i]-v2SD[Nend])/abs(v2SD[1]-v2SD[Nend]),
          abs(u3M[i]-u3M[Nend])/abs(u3M[1]-u3M[Nend]), abs(u3SD[i]-u3SD[Nend])/abs(u3SD[1]-u3SD[Nend]),
          abs(v3M[i]-v3M[Nend])/abs(v3M[1]-v3M[Nend]), abs(v3SD[i]-v3SD[Nend])/abs(v3SD[1]-v3SD[Nend])
  }
}
