 %clc;close all; 
load('codebook_J8_K4_200.mat','C')    % Load the codebook
M=size(C,2);  % size of the codebook for one user
m=log2(M); % no of bits in a symbol/block
K=4;   % No of orthogonal PREs/ OFDM subcarriers
J=8;   % No of users/layers
Nr=16; % No of BS antenaas
Nt=4;
ma=log2(Nt);
F=get_indicator_matrix_factor_graph(C, J, K);
d_v=length(find(F(:,1)));
d_f=length(find(F(1,:)));
%% %%%%%%% Power Allocation %%%%%%%%%
power_users=ones(1,J);  
sqrt_power_matrix=ones(K,J);
%% Monte Carlo %%%%%%%%%%%%%%%%%%%%%%%%
Eb_N0_dB=18;
%Eb_N0_dB=5:3:17;
%Eb_N0_dB=-5:2:5;
%Eb_N0_dB=18;
Eb_N0=10.^(Eb_N0_dB/10);
Es=sum(sum((abs(C)).^2))/length(C);
Eb=Es/m;   
N0=Eb./Eb_N0;
SNR=Es./N0;
SNR_db=10*log10(SNR);
sigma=sqrt(N0/2);
max_block_errors=[100 100 100 100 100 100 100 100 100 100 100];
%max_block_errors=[300  230 180 100 80 50 ];
%max_block_errors= 100   ;
SER=zeros(1,length(Eb_N0_dB));
%al=[0.3, 0.2, 0.4,0.3,0.4,0.7, 0.4,0.6,0.7,0.6];
al=[10, 25, 29,20,25];
T=15;
r=0:10:100;
%% 
 %for e=1:length(Eb_N0_dB)
 for e=1:length(r)
    rho=r(e)*N0/Es;
    block_errors=0;
    symbol_errors=0;
    block=0;
    while block_errors<max_block_errors(e)
        %%   SCMA Encoding %%%%%%

        bits_A1=randi([0 1],J,m);% blocks of bits for all users Tx-ant-1
        bits_A2=randi([0 1],J,m);% blocks of bits for all users Tx-ant-2
        bits_A3=randi([0 1],J,m);% blocks of bits for all users Tx-ant-3
        bits_A4=randi([0 1],J,m);% blocks of bits for all users Tx-ant-4
%         bits_A5=randi([0 1],J,m);% blocks of bits for all users Tx-ant-1
%         bits_A6=randi([0 1],J,m);% blocks of bits for all users Tx-ant-2
%         bits_A7=randi([0 1],J,m);% blocks of bits for all users Tx-ant-3
%         bits_A8=randi([0 1],J,m);% blocks of bits for all users Tx-ant-4
        %mapping
        symbol_indices_A1=bi2de(bits_A1,'left-msb')+1; % symbols for all users Tx-1
        symbol_indices_A2=bi2de(bits_A2,'left-msb')+1; % symbols for all users Tx-2
        symbol_indices_A3=bi2de(bits_A3,'left-msb')+1; % symbols for all users Tx-3
        symbol_indices_A4=bi2de(bits_A4,'left-msb')+1; % symbols for all users Tx-2
%         
%         symbol_indices_A5=bi2de(bits_A5,'left-msb')+1; % symbols for all users Tx-1
%         symbol_indices_A6=bi2de(bits_A6,'left-msb')+1; % symbols for all users Tx-2
%         symbol_indices_A7=bi2de(bits_A7,'left-msb')+1; % symbols for all users Tx-3
%         symbol_indices_A8=bi2de(bits_A8,'left-msb')+1; % symbols for all users Tx-2
        
        SCMA_codewords_A1=zeros(K,J);  % collection of the codewords for all users Tx-1 
       SCMA_codewords_A2=zeros(K,J);  % collection of the codewords for all users Tx-2
        SCMA_codewords_A3=zeros(K,J);  % collection of the codewords for all users Tx-3
        SCMA_codewords_A4=zeros(K,J);  % collection of the codewords for all users Tx-4 
%         
%         SCMA_codewords_A5=zeros(K,J);  % collection of the codewords for all users Tx-1 
%         SCMA_codewords_A6=zeros(K,J);  % collection of the codewords for all users Tx-2
%         SCMA_codewords_A7=zeros(K,J);  % collection of the codewords for all users Tx-3
%         SCMA_codewords_A8=zeros(K,J); 
        al_R=zeros(1,J);
        be_I=zeros(1,J);
        for j=1:J         % for each user
            present_codebook=C((j-1)*K+1:j*K,:);   % codebook for the jth user
            SCMA_codewords_A1(:,j)=present_codebook(:,symbol_indices_A1(j));
            SCMA_codewords_A2(:,j)=present_codebook(:,symbol_indices_A2(j));
            SCMA_codewords_A3(:,j)=present_codebook(:,symbol_indices_A3(j));
            SCMA_codewords_A4(:,j)=present_codebook(:,symbol_indices_A4(j));
%             
%             SCMA_codewords_A5(:,j)=present_codebook(:,symbol_indices_A5(j));
%             SCMA_codewords_A6(:,j)=present_codebook(:,symbol_indices_A6(j));
%             SCMA_codewords_A7(:,j)=present_codebook(:,symbol_indices_A7(j));
%             SCMA_codewords_A8(:,j)=present_codebook(:,symbol_indices_A8(j));
            al_R(1,j)=max(max(abs(real(present_codebook))));
            be_I(1,j)=max(max(abs(imag(present_codebook))));
%             
        end
        symbol_indices=[symbol_indices_A1 symbol_indices_A2 symbol_indices_A3 symbol_indices_A4];
%          symbol_indices=[symbol_indices_A1 symbol_indices_A2...
%                         symbol_indices_A3 symbol_indices_A4...
%                         symbol_indices_A5 symbol_indices_A6...
%                         symbol_indices_A7 symbol_indices_A8];
        %symbol_indices=[symbol_indices_A1 symbol_indices_A2 ];
           %symbol_indices=symbol_indices_A1;
               
        %% Transmission through Rayleigh fading channel %%
        AWGN=sigma*(randn(Nr,1)+1j*randn(Nr,1));  
        % complex Gaussian noise
        %AWGN=0;
        H_1=1/sqrt(2)*(randn(Nr,Nt*d_f)+1j*randn(Nr,Nt*d_f)); %RE-1
        H_2=1/sqrt(2)*(randn(Nr,Nt*d_f)+1j*randn(Nr,Nt*d_f)); %RE-2
        H_3=1/sqrt(2)*(randn(Nr,Nt*d_f)+1j*randn(Nr,Nt*d_f));%RE-3
        H_4=1/sqrt(2)*(randn(Nr,Nt*d_f)+1j*randn(Nr,Nt*d_f)); %RE-4
        %H_5=1/sqrt(2)*(randn(Nr,Nt*d_f)+1j*randn(Nr,Nt*d_f)); %RE-4 200%
        
        
    
         RE_CW=zeros(Nt*d_f,K);
        for k=1:K
            a=find(F(k,:));
            for df=1:d_f
              pr_RE_CW=[SCMA_codewords_A1(k,a(df)); SCMA_codewords_A2(k,a(df)); ...
                         SCMA_codewords_A3(k,a(df)); SCMA_codewords_A4(k,a(df))];                   
%                         SCMA_codewords_A5(k,a(df)) ;SCMA_codewords_A6(k,a(df));...
%                         SCMA_codewords_A7(k,a(df)) ;SCMA_codewords_A8(k,a(df))] ;
          % pr_RE_CW=[SCMA_codewords_A1(k,a(df)) ;SCMA_codewords_A2(k,a(df))];
           % pr_RE_CW=SCMA_codewords_A1(k,a(df)) ;
            RE_CW((df-1)*Nt+1:df*Nt,k)=pr_RE_CW;
            end
        end
        y_1=H_1*RE_CW(:,1)+AWGN;
        y_2=H_2*RE_CW(:,2)+AWGN;
        y_3=H_3*RE_CW(:,3)+AWGN;
        y_4=H_4*RE_CW(:,4)+AWGN; 
       % y_5=H_5*RE_CW(:,5)+AWGN;
        
           
     
        %% received SCMA codeword UP_link               
%          SCMA_CW_RE1_det=DPS_ADMM_SCMA_SMX_MIMO_150(H_1,y_1,rho,al,T,d_f,Nt,al_R,be_I);
%          SCMA_CW_RE2_det=DPS_ADMM_SCMA_SMX_MIMO_150(H_2,y_2,rho,al,T,d_f,Nt,al_R,be_I);
%          SCMA_CW_RE3_det=DPS_ADMM_SCMA_SMX_MIMO_150(H_3,y_3,rho,al,T,d_f,Nt,al_R,be_I);
%          SCMA_CW_RE4_det=DPS_ADMM_SCMA_SMX_MIMO_150(H_4,y_4,rho,al,T,d_f,Nt,al_R,be_I);
       a1=find(F(1,:));  
       al_R1=al_R(a1);be_I1=be_I(a1);
       SCMA_CW_RE1_det=DPS_ADMM_SCMA_SMX_MIMO_150(H_1,y_1,rho,al,T,d_f,Nt,al_R1,be_I1);
       a2=find(F(2,:));
       al_R2=al_R(a2);be_I2=be_I(a2);
       SCMA_CW_RE2_det=DPS_ADMM_SCMA_SMX_MIMO_150(H_2,y_2,rho,al,T,d_f,Nt,al_R2,be_I2);
       a3=find(F(3,:));
       al_R3=al_R(a3);be_I3=be_I(a3);
       SCMA_CW_RE3_det=DPS_ADMM_SCMA_SMX_MIMO_150(H_3,y_3,rho,al,T,d_f,Nt,al_R3,be_I3);
       a4=find(F(4,:));
       al_R4=al_R(a4);be_I4=be_I(a4);
       SCMA_CW_RE4_det=DPS_ADMM_SCMA_SMX_MIMO_150(H_4,y_4,rho,al,T,d_f,Nt,al_R4,be_I4);
%         a5=find(F(5,:));
%         al_R5=al_R(a5);be_I5=be_I(a5);
%         SCMA_CW_RE5_det=DPS_ADMM_SCMA_SMX_MIMO_150(H_5,y_5,rho,al,T,d_f,Nt,al_R5,be_I5);
%          
         
         
         
         
         
         
        %SCMA_CW_RE5_det=DPS_ADMM_SCMA_SM_150(H_5,y_5,rho,al,T,Nt,d_f);%200% overloading
%          SCMA_CW_RE1_det=MMSE_SCMA_SMX_150(H_1,y_1,Nt,d_f,rho);
%          SCMA_CW_RE2_det=MMSE_SCMA_SMX_150(H_2,y_2,Nt,d_f,rho);
%          SCMA_CW_RE3_det=MMSE_SCMA_SMX_150(H_3,y_3,Nt,d_f,rho);
%          SCMA_CW_RE4_det=MMSE_SCMA_SMX_150(H_4,y_4,Nt,d_f,rho);
        SCMA_CW_det=[ SCMA_CW_RE1_det SCMA_CW_RE2_det SCMA_CW_RE3_det SCMA_CW_RE4_det   ];
        
        SCMA_CW_det_mod=[reshape(SCMA_CW_RE1_det,[Nt,d_f]) reshape(SCMA_CW_RE2_det,[Nt,d_f])...
                         reshape(SCMA_CW_RE3_det,[Nt,d_f]) reshape(SCMA_CW_RE4_det,[Nt,d_f])];
               
        
       SCMA_Codewords_det_A1=zeros(K,J);
        for k=1:K
            a1=find(F(k,:));
           SCMA_Codewords_det_A1(k,[a1(1) a1(2) a1(3) a1(4)])=SCMA_CW_det_mod(1,(k-1)*d_f+1:k*d_f);
           %SCMA_Codewords_det_A1(k,[a1(1) a1(2) a1(3)])=SCMA_CW_det_mod(1,(k-1)*d_f+1:k*d_f);
        end  
        SCMA_Codewords_det_A2=zeros(K,J);
        for k=1:K
            a1=find(F(k,:));
           SCMA_Codewords_det_A2(k,[a1(1) a1(2) a1(3) a1(4)])=SCMA_CW_det_mod(2,(k-1)*d_f+1:k*d_f);
           %SCMA_Codewords_det_A2(k,[a1(1) a1(2) a1(3)])=SCMA_CW_det_mod(2,(k-1)*d_f+1:k*d_f);
        end  
        SCMA_Codewords_det_A3=zeros(K,J);
        for k=1:K
            a1=find(F(k,:));
           SCMA_Codewords_det_A3(k,[a1(1) a1(2) a1(3) a1(4)])=SCMA_CW_det_mod(3,(k-1)*d_f+1:k*d_f);
           %SCMA_Codewords_det_A3(k,[a1(1) a1(2) a1(3)])=SCMA_CW_det_mod(3,(k-1)*d_f+1:k*d_f);
        end  
        SCMA_Codewords_det_A4=zeros(K,J);
        for k=1:K
            a1=find(F(k,:));
           SCMA_Codewords_det_A4(k,[a1(1) a1(2) a1(3) a1(4)])=SCMA_CW_det_mod(4,(k-1)*d_f+1:k*d_f);
           %SCMA_Codewords_det_A4(k,[a1(1) a1(2) a1(3) ])=SCMA_CW_det_mod(4,(k-1)*d_f+1:k*d_f);
        end 
%           SCMA_Codewords_det_A5=zeros(K,J);
%         for k=1:K
%             a1=find(F(k,:));
%           SCMA_Codewords_det_A5(k,[a1(1) a1(2) a1(3)])=SCMA_CW_det_mod(5,(k-1)*d_f+1:k*d_f);
%         end  
%         SCMA_Codewords_det_A6=zeros(K,J);
%         for k=1:K
%             a1=find(F(k,:));
%            SCMA_Codewords_det_A6(k,[a1(1) a1(2) a1(3)])=SCMA_CW_det_mod(6,(k-1)*d_f+1:k*d_f);
%         end  
%         SCMA_Codewords_det_A7=zeros(K,J);
%         for k=1:K
%             a1=find(F(k,:));
%            SCMA_Codewords_det_A7(k,[a1(1) a1(2) a1(3)])=SCMA_CW_det_mod(7,(k-1)*d_f+1:k*d_f);
%         end  
%         SCMA_Codewords_det_A8=zeros(K,J);
%         for k=1:K
%             a1=find(F(k,:));
%            SCMA_Codewords_det_A8(k,[a1(1) a1(2) a1(3)])=SCMA_CW_det_mod(8,(k-1)*d_f+1:k*d_f);
%         end  
% % %         
% %         
%         
%         
%         
%         
%           
%         
%         
        det_ind_A1=zeros(J,1);
        for j=1:J           
             present_codebook=C((j-1)*K+1:j*K,:); %j-th user codebook
             ED1=zeros(1,M);
             for m1=1:M
                 ED=norm(SCMA_Codewords_det_A1(:,j)- present_codebook(:,m1));
                 ED1(1,m1)=ED;
             end
           det_ind_A1(j,1)=find(ED1==min(ED1));       
            %det_ind_A1(j,1)= det_ind_A1(j,1)-1;
            %det_bits_A1(j,:)=de2bi(det_ind_A1(j,1),'left-msb');
        end
        
        
        
        
        det_ind_A2=zeros(J,1);
        for j=1:J           
             present_codebook=C((j-1)*K+1:j*K,:); %j-th user codebook
             ED2=zeros(1,M);
             for m1=1:M
                 ED=norm(SCMA_Codewords_det_A2(:,j)- present_codebook(:,m1));
                 ED2(1,m1)=ED;
             end
           det_ind_A2(j,1)=find(ED2==min(ED2));       
              
        end
         det_ind_A3=zeros(J,1);
        for j=1:J           
             present_codebook=C((j-1)*K+1:j*K,:); %j-th user codebook
             ED3=zeros(1,M);
             for m1=1:M
                 ED=norm(SCMA_Codewords_det_A3(:,j)- present_codebook(:,m1));
                 ED3(1,m1)=ED;
             end
           det_ind_A3(j,1)=find(ED3==min(ED3));       
              
        end
        det_ind_A4=zeros(J,1);
        for j=1:J           
             present_codebook=C((j-1)*K+1:j*K,:); %j-th user codebook
             ED4=zeros(1,M);
             for m1=1:M
                 ED=norm(SCMA_Codewords_det_A4(:,j)- present_codebook(:,m1));
                 ED4(1,m1)=ED;
             end
           det_ind_A4(j,1)=find(ED4==min(ED4));       
              
        end
%         
%         det_ind_A5=zeros(J,1);
%         for j=1:J           
%              present_codebook=C((j-1)*K+1:j*K,:); %j-th user codebook
%              ED5=zeros(1,M);
%              for m1=1:M
%                  ED=norm(SCMA_Codewords_det_A5(:,j)- present_codebook(:,m1));
%                  ED5(1,m1)=ED;
%              end
%            det_ind_A5(j,1)=find(ED5==min(ED5));       
%               
%         end
%         
%           det_ind_A6=zeros(J,1);
%         for j=1:J           
%              present_codebook=C((j-1)*K+1:j*K,:); %j-th user codebook
%              ED6=zeros(1,M);
%              for m1=1:M
%                  ED=norm(SCMA_Codewords_det_A6(:,j)- present_codebook(:,m1));
%                  ED6(1,m1)=ED;
%              end
%            det_ind_A6(j,1)=find(ED6==min(ED6));       
%               
%         end
%         
%         det_ind_A7=zeros(J,1);
%         for j=1:J           
%              present_codebook=C((j-1)*K+1:j*K,:); %j-th user codebook
%              ED7=zeros(1,M);
%              for m1=1:M
%                  ED=norm(SCMA_Codewords_det_A7(:,j)- present_codebook(:,m1));
%                  ED7(1,m1)=ED;
%              end
%            det_ind_A7(j,1)=find(ED7==min(ED7));       
%               
%         end
%         
%         det_ind_A8=zeros(J,1);
%         for j=1:J           
%              present_codebook=C((j-1)*K+1:j*K,:); %j-th user codebook
%              ED8=zeros(1,M);
%              for m1=1:M
%                  ED=norm(SCMA_Codewords_det_A8(:,j)- present_codebook(:,m1));
%                  ED8(1,m1)=ED;
%              end
%            det_ind_A8(j,1)=find(ED8==min(ED8));       
%               
%         end
        
        
        
        
        
        %det_ind_all_A=det_ind_A1;
        det_ind_all_A=[det_ind_A1 det_ind_A2 det_ind_A3 det_ind_A4  ];
        %det_ind_all_A=[det_ind_A1 det_ind_A2 det_ind_A3 det_ind_A4 det_ind_A5 det_ind_A6 det_ind_A7 det_ind_A8  ];
        %number_errors=sum(or(det_ind_A1~=symbol_indices_A1,det_ind_A2~= symbol_indices_A2,...
            %det_ind_A3~= symbol_indices_A3,det_ind_A4~= symbol_indices_A4));
        %error_locations=find(det_ind~=symbol_indices);
        %if  ~isempty(error_locations)
        errors_locations=find(symbol_indices~=det_ind_all_A);
        if ~isempty( errors_locations)
            block_errors=block_errors+1;           
            %symbol_errors=symbol_errors+length(error_locations);
            symbol_errors=symbol_errors+length( errors_locations);
           % fprintf('\n  %d error collected',block_errors); 
        end     
        block=block+1;
    end
    SER(e)=symbol_errors/block/J/Nt;% Each block contains J-symbols
    %BER(e)=SER(e)/m;
    fprintf('\n Simulation done for %d ',r(e)); 
end
semilogy(r,SER,'b-*','LineWidth',2) ;
xlabel('r');
ylabel('SER');
legend(' Nr=16, Nt=2 M=4 rho vs SER');





