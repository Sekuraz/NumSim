BEGIN {
  FS="\t";
  OFS="\t";
  PI=atan2(1,1)*4;
  MU=1500;
  SIGMA=166.666666666667;
  FACTOR=1/(SIGMA * sqrt(2*PI));
  FACTOR2=-1/(2*SIGMA^2);
  N=200;
}
function abs(x) { return (x<0) ? -x : x }
function gausspdf(x) {
  return FACTOR * exp(FACTOR2*(x-MU)^2);
}
{
  p=gausspdf($1);
  
  for(i=1; i<=N; i++) {
    if(N%i == 0 && (NR-1)%i == 0) {
      nr[i]++;
      psum[i]+=p;
      re[i]+=$1 * p; reSq[i]+=$1^2 * p;
      u1[i]+=$4 * p; u1Sq[i]+=$4^2 * p;
      v1[i]+=$5 * p; v1Sq[i]+=$5^2 * p;
      u2[i]+=$6 * p; u2Sq[i]+=$6^2 * p;
      v2[i]+=$7 * p; v2Sq[i]+=$7^2 * p;
      u3[i]+=$8 * p; u3Sq[i]+=$8^2 * p;
      v3[i]+=$9 * p; v3Sq[i]+=$9^2 * p;
    }
  }
}
END {
  for(i=1; i<=N; i++) {
    if(N%i == 0 && (NR-1)%i == 0) {
      reSD[i]=sqrt((reSq[i] - re[i]^2/psum[i])/psum[i]);
      u1SD[i]=sqrt((u1Sq[i] - u1[i]^2/psum[i])/psum[i]);
      v1SD[i]=sqrt((v1Sq[i] - v1[i]^2/psum[i])/psum[i]);
      u2SD[i]=sqrt((u2Sq[i] - u2[i]^2/psum[i])/psum[i]);
      v2SD[i]=sqrt((v2Sq[i] - v2[i]^2/psum[i])/psum[i]);
      u3SD[i]=sqrt((u3Sq[i] - u3[i]^2/psum[i])/psum[i]);
      v3SD[i]=sqrt((v3Sq[i] - v3[i]^2/psum[i])/psum[i]);
      reM[i]=re[i]/psum[i];
      u1M[i]=u1[i]/psum[i]; v1M[i]=v1[i]/psum[i];
      u2M[i]=u2[i]/psum[i]; v2M[i]=v2[i]/psum[i];
      u3M[i]=u3[i]/psum[i]; v3M[i]=v3[i]/psum[i];
    }
  }
  for(i=2; i<=N; i++) {
    if(N%i == 0 && (NR-1)%i == 0) {
      print nr[i],
          abs(reM[i]-1500),#/abs(reM[N]-1500),
          abs(reSD[i]-166.666666666667)/abs(reSD[N]-166.666666666667),
          abs(u1M[i]-u1M[1])/abs(u1M[N]-u1M[1]),
          abs(u1SD[i]-u1SD[1])/abs(u1SD[N]-u1SD[1]),
          abs(v1M[i]-v1M[1])/abs(v1M[N]-v1M[1]),
          abs(v1SD[i]-v1SD[1])/abs(v1SD[N]-v1SD[1]),
          abs(u2M[i]-u2M[1])/abs(u2M[N]-u2M[1]),
          abs(u2SD[i]-u2SD[1])/abs(u2SD[N]-u2SD[1]),
          abs(v2M[i]-v2M[1])/abs(v2M[N]-v2M[1]),
          abs(v2SD[i]-v2SD[1])/abs(v2SD[N]-v2SD[1]),
          abs(u3M[i]-u3M[1])/abs(u3M[N]-u3M[1]),
          abs(u3SD[i]-u3SD[1])/abs(u3SD[N]-u3SD[1]),
          abs(v3M[i]-v3M[1])/abs(v3M[N]-v3M[1]),
          abs(v3SD[i]-v3SD[1])/abs(v3SD[N]-v3SD[1])
    }
  }
}
