BEGIN {
  FS="\t";
  OFS="\t";
  PI=atan2(1,1)*4;
  MU=1500;
  SIGMA=166.666666666667;
  FACTOR=1/(SIGMA * sqrt(2*PI));
  FACTOR2=-1/(2*SIGMA^2);
  N=50;
}
function gausspdf(x) {
  return FACTOR * exp(FACTOR2*(x-MU)^2);
}
(FNR-1)%N == 0 {
  nr[FNR]++;
  p=gausspdf($1);
  psum[FNR]+=p;
  re[FNR]+=$1 * p; reSq[FNR]+=$1^2 * p;
  t[FNR]+=$2;
  res[FNR]+=$3 * p; resSq[FNR]+=$3^2 * p;
  u1[FNR]+=$4 * p; u1Sq[FNR]+=$4^2 * p;
  v1[FNR]+=$5 * p; v1Sq[FNR]+=$5^2 * p;
  u2[FNR]+=$6 * p; u2Sq[FNR]+=$6^2 * p;
  v2[FNR]+=$7 * p; v2Sq[FNR]+=$7^2 * p;
  u3[FNR]+=$8 * p; u3Sq[FNR]+=$8^2 * p;
  v3[FNR]+=$9 * p; v3Sq[FNR]+=$9^2 * p;
}
END {
  for(i=1; i<=FNR; i+=N) {
    reSD=sqrt((reSq[i] - (re[i])^2/psum[i])/psum[i]);
    resSD=sqrt((resSq[i] - (res[i])^2/psum[i])/psum[i]);
    u1SD=sqrt((u1Sq[i] - (u1[i])^2/psum[i])/psum[i]);
    v1SD=sqrt((v1Sq[i] - (v1[i])^2/psum[i])/psum[i]);
    u2SD=sqrt((u2Sq[i] - (u2[i])^2/psum[i])/psum[i]);
    v2SD=sqrt((v2Sq[i] - (v2[i])^2/psum[i])/psum[i]);
    u3SD=sqrt((u3Sq[i] - (u3[i])^2/psum[i])/psum[i]);
    v3SD=sqrt((v3Sq[i] - (v3[i])^2/psum[i])/psum[i]);
    print t[i]/nr[i],
        re[i]/psum[i], reSD,
        res[i]/psum[i], resSD,
        u1[i]/psum[i], u1SD, u1[i]/psum[i] + u1SD, u1[i]/psum[i] - u1SD,
        v1[i]/psum[i], v1SD, v1[i]/psum[i] + v1SD, v1[i]/psum[i] - v1SD,
        u2[i]/psum[i], u2SD, u2[i]/psum[i] + u2SD, u2[i]/psum[i] - u2SD,
        v2[i]/psum[i], v2SD, v2[i]/psum[i] + v2SD, v2[i]/psum[i] - v2SD,
        u3[i]/psum[i], u3SD, u3[i]/psum[i] + u3SD, u3[i]/psum[i] - u3SD,
        v3[i]/psum[i], v3SD, v3[i]/psum[i] + v3SD, v3[i]/psum[i] - v3SD
  }
}
