function SCMA_CW_RE_det=DPS_ADMM_SCMA_SMX_MIMO_150(H,y,rho,al,T,d_f,Nt,al_R,be_I)
I=eye(Nt*d_f)/((H'*H)*d_f+rho*eye(Nt*d_f));
Ma=H'*y;
%Initializing
x_11=zeros(Nt*d_f,1)+zeros(Nt*d_f,1)*1i;
x_12=zeros(Nt*d_f,1)+zeros(Nt*d_f,1)*1i;
x_13=zeros(Nt*d_f,1)+zeros(Nt*d_f,1)*1i;
x_14=zeros(Nt*d_f,1)+zeros(Nt*d_f,1)*1i;
 %x_15=zeros(d_f,1)+zeros(d_f,1)*1i;
% x_16=zeros(d_f,1)+zeros(d_f,1)*1i;
x_1_bar=zeros(Nt*d_f,1)+zeros(Nt*d_f,1)*1i;
x_0_bar=zeros(Nt*d_f,1)+zeros(Nt*d_f,1)*1i;
u=zeros(Nt*d_f,1)+zeros(Nt*d_f,1)*1i;
for t=1:T
   x_11=(rho/(rho+al(1)))*(x_11+x_0_bar-u-x_1_bar);
   x_12=(rho/(rho+al(2)))*(x_12+x_0_bar-u-x_1_bar);
   x_13=(rho/(rho+al(3)))*(x_13+x_0_bar-u-x_1_bar);
   x_14=(rho/(rho+al(4)))*(x_14+x_0_bar-u-x_1_bar);
%    x_15=(rho/(rho+al(5)))*(x_15+x_0_bar-u-x_1_bar);
%    x_16=(rho/(rho+al(6)))*(x_16+x_0_bar-u-x_1_bar);
   
   x_11R=real(x_11);x_11I=imag(x_11);
   x_11R(x_11R<-al_R(1))=-al_R(1); x_11R(x_11R>al_R(1))=al_R(1);
   x_11I(x_11I<-be_I(1))=-be_I(1);x_11I(x_11I>be_I(1))=be_I(1);
   x_12R=real(x_12);x_12I=imag(x_12);
   x_12R(x_12R<-al_R(2))=-al_R(2); x_12R(x_12R>al_R(2))=al_R(2);
   x_12I(x_12I<-be_I(2))=-be_I(2);x_12I(x_12I>be_I(2))=be_I(2);
   x_13R=real(x_13);x_13I=imag(x_13);
 x_13R(x_13R<-al_R(3))=-al_R(3); x_13R(x_13R>al_R(3))=al_R(3);
   x_13I(x_13I<-be_I(3))=-be_I(3);x_13I(x_13I>be_I(3))=be_I(3);
     x_14R=real(x_14);x_14I=imag(x_14);
   x_14R(x_14R<-al_R(4))=-al_R(4); x_14R(x_14R>al_R(4))=al_R(4);
%    x_14I(x_14I<-be_I(4))=-be_I(4);x_14I(x_14I>be_I(4))=be_I(4);
%    x_14R=real(x_14);x_14I=imag(x_14);
%    x_14R(x_14R<-1)=-1; x_14R(x_14R>1)=1;
%    x_14I(x_14I<-1)=-1;x_14I(x_14I>1)=1;
%    x_15R=real(x_15);x_15I=imag(x_15);
%    x_15R(x_15R<-1)=-1; x_15R(x_15R>1)=1;
%    x_15I(x_15I<-1)=-1;x_15I(x_15I>1)=1;
%    x_16R=real(x_16);x_16I=imag(x_16);
%    x_16R(x_16R<-1)=-1; x_16R(x_16R>1)=1;
%    x_16I(x_16I<-1)=-1;x_16I(x_16I>1)=1;  
   x_11=x_11R+x_11I*1j;x_12=x_12R+x_12I*1j;x_13=x_13R+x_13I*1j;
   x_14=x_14R+x_14I*1j;
   %x_15=x_15R+x_15I*1j;x_16=x_16R+x_16I*1j;  
   x_11(Nt+1:end)=0;x_12(1:Nt)=0;x_12(2*Nt+1:end)=0;x_13(1:2*Nt)=0;
   x_13(3*Nt+1:end)=0;
   x_14(1:3*Nt)=0;
   %x_14(1:6)=0;
   %x_14(9:end)=0;x_15(1:8)=0;x_15(11:end)=0;x_16(1:10)=0;
   x_1_bar=sum(x_11+x_12+x_13+x_14)/d_f;
   %x_1_bar=sum(x_11+x_12+x_13)/d_f;
   x_0_bar=I*(Ma+rho*(x_1_bar+u));
   u=u+(x_1_bar-x_0_bar);      
end
SCMA_CW_RE_det=d_f*x_0_bar;



