function dy = Fun(t,y,ka,kda,C,kc,ksb,krb,kra,Gj,Gd,Gb,Gt,Gg,Gh,Gl,sigmaj,sigmad,sigmab,sigmat,sigmag,sigmah,sigmal,kcasp,kin,kout,wbec0,wjnk0,wdapk0,wbcl0,watg0,wmtor0,wbh0,wjnk_s,wdapk_s,wbec_d,wbcl_j,watg_m,wmtor_ca,wbh_s,BCL2T,IP3RT,BCLmit,CaT,CASP,Wcp_0,Wcp_Ca,W_ATG50,W_ATG5C,Ga,Gcp,sigmacp,sigmaA,wbh_c,wbh_a)
BECNT = y(1);Atphg  = y(2);S = y(3);JNK  = y(4); DAPK  = y(5); BCLP = y(6); mTOR  = y(7); Atg13  = y(8); Ca  = y(9);BH3 = y(10); BECNP = y(11); calpain = y(12); Atg5 = y(13); 
 
%Algebric Equations
BCL2U=BCL2T-BCLP;
BECNU=BECNT-BECNP;
LIGtotal=IP3RT+BECNU; 
LIGfree= max(0,(LIGtotal-BCL2U));
BECNF=BECNU*LIGfree/LIGtotal;
IP3RF=IP3RT*LIGfree/LIGtotal;
% 
if (BH3-BCLmit)>0
   CASP = 1;
end
 
%ODEs
dy(1)=-kcasp*CASP*BECNT;
dy(2)=ka*(BECNF+BECNP)*Atg13*Atg5-kda*Atphg;
dy(3)=ksb+C*kc-(krb+kra*Atphg)*S;
dy(4)=Gj*(1/(1+exp(-sigmaj*(wjnk0+wjnk_s*S)))-JNK);
dy(5)=Gd*(1/(1+exp(-sigmad*(wdapk0+wdapk_s*S)))-DAPK);
dy(6)=Gb*(BCL2T*1/(1+exp(-sigmab*(wbcl0+wbcl_j*JNK)))-BCLP);
dy(7)=Gt*(1/(1+exp(-sigmat*(wmtor0+wmtor_ca*Ca)))-mTOR);
dy(8)=Gg*(1/(1+exp(-sigmag*(watg0+watg_m*mTOR)))-Atg13);
dy(9)=kout*IP3RF*(CaT-Ca)-kin*Ca;
dy(10)=Gh*(1/(1+exp(-sigmah*(wbh0+wbh_s*S + wbh_c*Ca + wbh_a*(1-Atg5))))-BH3);
dy(11)=Gl*(BECNT*1/(1+exp(-sigmal*(wbec0+wbec_d*DAPK)))-BECNP);
dy(12)=Gcp*(1/(1+exp(-sigmacp*(Wcp_0+Wcp_Ca*Ca)))-calpain);
dy(13)= Ga*(1/(1+exp(-sigmaA*(W_ATG50 + W_ATG5C*calpain)))-Atg5);
dy=dy';
