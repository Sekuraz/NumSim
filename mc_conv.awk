BEGIN {
  FS="\t";
  OFS="\t";
}
{
  nr++;
  re+=$1; reSq+=$1^2;
  u1+=$4; u1Sq+=$4^2; v1+=$5; v1Sq+=$5^2;
  u2+=$6; u2Sq+=$6^2; v2+=$7; v2Sq+=$7^2;
  u3+=$8; u3Sq+=$8^2; v3+=$9; v3Sq+=$9^2;
  
  if(nr % 500 == 0) {
    reSD=sqrt((reSq - re^2/nr)/(nr-1));
    u1SD=sqrt((u1Sq - u1^2/nr)/(nr-1));
    v1SD=sqrt((v1Sq - v1^2/nr)/(nr-1));
    u2SD=sqrt((u2Sq - u2^2/nr)/(nr-1));
    v2SD=sqrt((v2Sq - v2^2/nr)/(nr-1));
    u3SD=sqrt((u3Sq - u3^2/nr)/(nr-1));
    v3SD=sqrt((v3Sq - v3^2/nr)/(nr-1));
    print nr,
          re/nr, reSD,
          u1/nr, u1SD, v1/nr, v1SD,
          u2/nr, u2SD, v2/nr, v2SD,
          u3/nr, u3SD, v3/nr, v3SD
  }
}
