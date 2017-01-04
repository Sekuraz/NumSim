BEGIN {
  FS="\t";
  OFS="\t";
  PI=atan2(1,1)*4;
  MU=1500;
  SIGMA=166.666666666667;
  FACTOR=1/(SIGMA * sqrt(2*PI));
  FACTOR2=-1/(2*SIGMA^2);
}
function gausspdf(x) {
  return FACTOR * exp(FACTOR2*(x-MU)^2);
}
{
  nr++;
  p=gausspdf($1);
  psum+=p;
  re+=$1 * p; reSq+=$1^2 * p;
  u1+=$4 * p; u1Sq+=$4^2 * p;
  v1+=$5 * p; v1Sq+=$5^2 * p;
  u2+=$6 * p; u2Sq+=$6^2 * p;
  v2+=$7 * p; v2Sq+=$7^2 * p;
  u3+=$8 * p; u3Sq+=$8^2 * p;
  v3+=$9 * p; v3Sq+=$9^2 * p;
  
  if(nr % 50 == 0) {
    reSD=sqrt((reSq - re^2/psum)/psum);
    u1SD=sqrt((u1Sq - u1^2/psum)/psum);
    v1SD=sqrt((v1Sq - v1^2/psum)/psum);
    u2SD=sqrt((u2Sq - u2^2/psum)/psum);
    v2SD=sqrt((v2Sq - v2^2/psum)/psum);
    u3SD=sqrt((u3Sq - u3^2/psum)/psum);
    v3SD=sqrt((v3Sq - v3^2/psum)/psum);
    print nr,
          re/psum, reSD,
          u1/psum, u1SD, v1/psum, v1SD,
          u2/psum, u2SD, v2/psum, v2SD,
          u3/psum, u3SD, v3/psum, v3SD
  }
}
