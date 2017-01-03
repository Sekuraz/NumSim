BEGIN {
  FS="\t";
  OFS="\t";
}
{
  nr[FNR]++;
  re[FNR]+=$1; reSq[FNR]+=$1^2;
  t[FNR]+=$2;
  res[FNR]+=$3; resSq[FNR]+=$3^2;
  u1[FNR]+=$4; u1Sq[FNR]+=$4^2; v1[FNR]+=$5; v1Sq[FNR]+=$5^2;
  u2[FNR]+=$6; u2Sq[FNR]+=$6^2; v2[FNR]+=$7; v2Sq[FNR]+=$7^2;
  u3[FNR]+=$8; u3Sq[FNR]+=$8^2; v3[FNR]+=$9; v3Sq[FNR]+=$9^2;
}
END {
  for(i=1; i<=FNR; i++) {
    reSD=sqrt((reSq[i] - (re[i])^2/nr[i])/(nr[i]-1));
    resSD=sqrt((resSq[i] - (res[i])^2/nr[i])/(nr[i]-1));
    u1SD=sqrt((u1Sq[i] - (u1[i])^2/nr[i])/(nr[i]-1));
    v1SD=sqrt((v1Sq[i] - (v1[i])^2/nr[i])/(nr[i]-1));
    u2SD=sqrt((u2Sq[i] - (u2[i])^2/nr[i])/(nr[i]-1));
    v2SD=sqrt((v2Sq[i] - (v2[i])^2/nr[i])/(nr[i]-1));
    u3SD=sqrt((u3Sq[i] - (u3[i])^2/nr[i])/(nr[i]-1));
    v3SD=sqrt((v3Sq[i] - (v3[i])^2/nr[i])/(nr[i]-1));
    print t[i]/nr[i],
        re[i]/nr[i], reSD,
        res[i]/nr[i], resSD,
        u1[i]/nr[i], u1SD, u1[i]/nr[i] + u1SD, u1[i]/nr[i] - u1SD,
        v1[i]/nr[i], v1SD, v1[i]/nr[i] + v1SD, v1[i]/nr[i] - v1SD,
        u2[i]/nr[i], u2SD, u2[i]/nr[i] + u2SD, u2[i]/nr[i] - u2SD,
        v2[i]/nr[i], v2SD, v2[i]/nr[i] + v2SD, v2[i]/nr[i] - v2SD,
        u3[i]/nr[i], u3SD, u3[i]/nr[i] + u3SD, u3[i]/nr[i] - u3SD,
        v3[i]/nr[i], v3SD, v3[i]/nr[i] + v3SD, v3[i]/nr[i] - v3SD
  }
}
