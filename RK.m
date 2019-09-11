clc;
clear all;

disp('========== Simulasi Sistem Konduksi Panas ==========')
disp('Tahapan 2 : Sistem Direduksi Kemudian Dikontrol')
disp(' ')

deltaT = input('Masukkan lama waktu perambatan (Satuan Sekon) = ');
deltaX = input ('Masukkan panjang perambatan panas (Satuan meter) = ');
fprintf('Jenis Konduktor : \n1. Emas \n2. Perak \n3. Tembaga \n4. Aluminium \n5. Baja \n6. Besi \n');
alfa = input('Pilih Jenis Konduktor = ');
%Temperatur satuan celcius
switch alfa
    case 1
        alfa = 300;
    case 2
        alfa = 420;
    case 3
        alfa = 380;
    case 4
        alfa = 200;
    case 5
        alfa = 40;
    case 6
        alfa = 80;
end 

a = alfa*deltaT/(deltaX)^2;
%a = 0.4;
n = input('Masukkan orde yang diinginkan = ');

%Kontruksi matriks A
A=zeros(n);
for i=1:n
    for j=1:n
        if j==i
            A(i,j)=1-2*a;
        elseif i-j==-1 | i-j==1 
            A(i,j)=a;
        end
    end 

end
%Kontriksi matriks B
B = zeros(n,1);
B(1,1)=a;
%Kontruksi matriks C
C = zeros(1,n);
C(1,1)=1;
%Kontruksi matriks D
D=[1];

tic
disp('========= Sistem Awal =========')
sys=ss(A,B,C,D,1)
% Fungsi transfer sistem awal
G_awal = tf(sys); 

%% Analsisa Sistem Awal Konduksi Panas

disp('=========================================')
disp(' ')
disp('Analisa Sistem Awal')
disp(' ')
disp('Nilai Eigen Dari Sistem :')
Eigen=abs(eig(A))
%inisiasi, jumlah eigen stabil dan tidak stabil
Tidak_Stabil = 0;
Stabil = 0;
Stabil_Asimtotik = 0;
n = size(A);
for i = 1:n
    if real(Eigen(i)) > 1
       Tidak_Stabil = Tidak_Stabil +1;
    end
        if real(Eigen(i)) < 1
            Stabil_Asimtotik = Stabil_Asimtotik +1;
        end
            if real(Eigen(i)) == 1
               Stabil = Stabil +1;
            end
end
Stabil=Stabil+Stabil_Asimtotik;
Tidak_Stabil;
Stabil;

if Stabil == n
    disp('Maka Sistem Stabil')
else
    disp('Maka Sistem Tidak Stabil')
end
disp(' ')
disp('=========================================')
fprintf('Uji Keterkontrolan dan Uji Keteramatan Sistem Awal \n\n')
%Keterkontrolan
Mc = ctrb(sys);
RankMc = rank(Mc);
%Ketermatan
Mo = obsv(sys);
RankMo = rank(Mo);

if rank(A)==rank(Mc)
    if rank(A)==rank(Mo)
        fprintf('rank A = rank Mc = rank Mo = %.0f\n',n');
        fprintf('Maka Sistem terkontrol dan teramati\n');
    else
        fprintf('rank A = rank Mc, rank A =/= rank Mo\n\n');
        fprintf('Maka Sistem terkontrol namun tidak teramati\n');
    end
else rank(A)~= rank(Mc)
    if rank(A)==rank(Mo)
        fprintf('rank A =/= rank Mc, rank A = rank Mo\n\n');
        fprintf('Maka Sistem tidak terkontrol namun teramati\n');
    else
        fprintf('rank A =/= rank Mc =/= rank Mo\n\n');
        fprintf('Maka Sistem tidak terkontrol dan tidak teramati\n');
    end
end

%'=========================================')
%Gramian Keterkontrolan
W = gram(sys,'c');
%Cek W definit positif
EigenW = eig(W);
%Cek W nonsingular
DeterminanW = det(W);
%'=========================================')
%Gramian Keteramatan
M = gram(sys,'o');
%Cek M definit positif
EigenW = eig(M);
%Cek M nonsingular
DeterminanW = det(M);
disp(' ')


%% Kontruksi matriks T

%Menentukan matriks phi dimana berlaku W=phi(transpose)*phi'
phi=chol(W);

%Cek W=phi(transpose)*phi
Cek = phi'*phi;

%Diagonalisai phi*M*phi(transpose) sedemikian hingga berlaku phi*M*phi(transpose) = U*(sigma^2)*U(transpose)
J = phi*M*phi';
[U,S,V]=svd(J);

%muncul J=U*S*V'
%S = sigma^2
%yang dicari diagonalisai J sehingga J=U*(sigma^2)*U'
%sehingga akar S menunjukkan nilai singular hankel
sigma = (S).^(1/2);

%Diperoleh Matriks Transformasi T
T=phi'*U*(sigma)^(-0.5);

%% Kontruksi Matriks Setimbang

disp('=========================================')
disp('Sistem Setimbang Dari Konduksi Panas')
disp(' ')
As=(inv(T))*A*T;
Bs=(inv(T))*B;
Cs=C*T;
Ds=D;
sys_s=ss(As,Bs,Cs,Ds,1)

%Fungsi transfer sistem setimbang
G_setimbang = tf(sys_s);


%% Analisa Sistem Setimbang Konduksi Panas

disp('=========================================')
disp('========= Analisa Sistem Setimbang =========')
disp(' ')
disp('Cek Kestabilan Sistem Setimbang')
Eigen_S=abs(eig(As))
%inisiasi, jumlah eigen stabil dan tidak stabil
Tidak_Stabil_S = 0;
Stabil_S = 0;
Stabil_Asimptotik_S = 0;
n = size(As);
for i = 1:n
    if real(Eigen_S(i)) > 1
       Tidak_Stabil_S = Tidak_Stabil_S +1;
    end
        if real(Eigen_S(i)) < 1
            Stabil_Asimptotik_S = Stabil_Asimptotik_S +1;
        end
            if real(Eigen_S(i)) == 1
               Stabil_S = Stabil_S +1;
            end
end
Stabil_S=Stabil_S+Stabil_Asimptotik_S;
Tidak_Stabil_S;
Stabil_S;

if Stabil_S == n
    disp('Maka Sistem Stabil')
else
    disp('Maka Sistem Tidak Stabil')
end

disp('')
disp('=========================================')
fprintf('Uji Keterkontrolan dan Uji Keteramatan Sistem Setimbang \n\n')
%Keterkontrolan
Mc_s = ctrb(sys_s);
RankMc_s = rank(Mc_s);
%Ketermatan
Mo_s = obsv(sys_s);
RankMo_s = rank(Mo_s);

if rank(A)==rank(Mc)
    if rank(A)==rank(Mo)
        fprintf('rank A = rank Mc = rank Mo = %.0f\n',n');
        fprintf('Maka Sistem, terkontrol dan teramati\n');
    else
        fprintf('rank A = rank Mc , rank A =/= rank Mo \n\n');
        fprintf('Maka Sistem Setimbang, terkontrol namun tidak teramati\n');
    end
else rank(A)~= rank(Mc)
    if rank(A)==rank(Mo)
        fprintf('rank A =/= rank Mc , rank A  = rank Mo  \n\n');
       fprintf('Maka Sistem Setimbang, tidak terkontrol namun teramati\n');
    else
       fprintf('rank A =/= rank Mc =/= rank Mo \n\n');
       fprintf('Maka Sistem Setimbang, tidak terkontrol dan tidak teramati\n');
    end
end

disp('')
disp('=========================================')
disp('Gramian Keterkontrolan Sistem Setimbang')
Ws = gram(sys_s,'c')

disp('=========================================')
disp('Gramian Keteramatan Sistem Setimbang')
Ms = gram(sys_s,'o')

disp('=========================================')
disp('Nilai Singular Hankel Sistem Setimbang')
hsv = hsvd(sys_s);
disp(hsv)


%% Reduksi Sistem Setimbang 

k = size(As);
disp(['Masukkan ukuran sistem tereduksi yang diinginkan (ukuran tidak lebih dari orde matriks A)']);
r = input('r = '); 

%partisi matriks sistem setimbang
A11 = As(1:r,1:r);
A12 = As(1:r,r+1:k);
A21 = As(r+1:k,1:r);
A22 = As(r+1:k,r+1:k);

B1 = Bs(1:r);
B2 = Bs(r+1:k);

C1 = Cs(1:r);
C2 = Cs(1+r:k);

Sigma1 = sigma(1:r,1:r);
Sigma2 = sigma(r+1,r+1);
sigma_r = 0;
for i=1+r:k
    sigma_k = hsv(i);
    sigma_j = 2*sigma_k + sigma_r;
    sigma_r = sigma_j;
end

Ar = A11;
Br = B1;
Cr = C1;
Dr = D;

Sys_Tereduksi = ss(Ar,Br,Cr,Dr,1)
G_TereduksiBT = tf(Sys_Tereduksi);

%% Analisa Sistem Tereduksi


disp('=========================================')
disp('========= Analisa Sistem Tereduksi =========')
disp(' ')
disp('Nilai Eigen Sistem Tereduksi')
Eigen_tereduksi = abs(eig(Ar))
%inisiasi, jumlah eigen stabil dan tidak stabil
Tidak_Stabil_R = 0;
Stabil_R = 0;
Stabil_Asimptotik_R = 0;
n = size(Ar);
for i = 1:n
    if real(Eigen_tereduksi(i)) > 1
       Tidak_Stabil_R = Tidak_Stabil_R +1;
    end
        if real(Eigen_tereduksi(i)) < 1
            Stabil_Asimptotik_R = Stabil_Asimptotik_R +1;
        end
            if real(Eigen_tereduksi(i)) == 1
               Stabil_R = Stabil_R +1;
            end
end
Stabil_R=Stabil_R+Stabil_Asimptotik_R;
Tidak_Stabil_R;
Stabil_R;

if Stabil_R == n
    disp('Maka Sistem Stabil')
else
    disp('Maka Sistem Tidak Stabil')
end

disp('')
disp('=========================================')
fprintf('Uji Keterkontrolan dan Uji Keteramatan Sistem Tereduksi \n\n')
%Keterkontrolan
Mc_Tereduksi = ctrb(Sys_Tereduksi);
RankMc_Tereduksi = rank(Mc_Tereduksi);
%Ketermatan
Mo_Tereduksi = obsv(Sys_Tereduksi);
RankMo_Tereduksi = rank(Mo_Tereduksi);

if rank(A)==rank(Mc)
    if rank(A)==rank(Mo)
        fprintf('rank A = rank Mc = rank Mo = %.0f\n',n');
        fprintf('Maka Sistem Tereduksi, terkontrol dan teramati\n');
    else
        fprintf('rank A = rank Mc , rank A =/= rank Mo  \n\n');
        fprintf('Maka Sistem Tereduksi, terkontrol namun tidak teramati\n');
    end
else rank(A)~= rank(Mc)
    if rank(A)==rank(Mo)
       fprintf('rank A =/= rank Mc , rank A = rank Mo \n\n');
       fprintf('Maka Sistem Tereduksi, tidak terkontrol namun teramati\n');
    else
       fprintf('rank A =/= rank Mc =/= rank Mo \n\n');
       fprintf('Maka Sistem Tereduksi, tidak terkontrol dan tidak teramati\n');
    end
end

disp('')

%Gramian Keterkontrolan
W_Tereduksi = gram(Sys_Tereduksi,'c');
%Cek W definit positif
EigenW_Tereduksi = eig(W_Tereduksi);
%Cek W nonsingular
DeterminanW_Tereduksi = det(W_Tereduksi);

%Gramian Keteramatan
M_Tereduksi = gram(Sys_Tereduksi,'o');
%Cek M definit positif
EigenW_Tereduksi = eig(M_Tereduksi);
%Cek M nonsingular
DeterminanW_Tereduksi = det(M_Tereduksi);

disp('=========================================')
%Cek eror sistem tereduksi
Error_BT=G_awal-G_TereduksiBT;

disp('Cek Syarat Reduksi Dengan Balanced Truncation Norm_BT =< NilaiSigma_r')
Norm_BT=norm(Error_BT,inf)
NilaiSigma_r = sigma_r
if Norm_BT <= NilaiSigma_r
    disp('Syarat Terpenuhi')
else
    disp('Syarat Tidak Terpenuhi')
end
disp('')
disp('Nilai error sistem tereduksi dengan Balanced Truncation : ')
disp(Norm_BT)

%% Desain Kontrol Dengan LQR
%Kontrol Sistem Tereduksi

disp(' ')
disp('======== Sistem Tereduksi Yang Dikontrol ========')
Qs = Cr'*Cr
Rs=1;
[Ks, Ss, Es] = dlqr(Ar,Br,Qs,Rs)%Ss=Penyelesaian Pers. Riccati, Es = Eigenvalue
Ar_C=Ar-Br*Ks;
Br_C=Br;
Cr_C=Cr;
Dr_C=Dr;
sys_reduksi_kontrol=ss(Ar_C,Br_C,Cr_C,Dr_C,1)
G_Tereduksi_Kontrol = tf(sys_reduksi_kontrol);


%% Analsisa Sistem Tereduksi Yang Dikontrol

disp('=========================================')
disp('Analisa Sistem Tereduksi Yang Dikontrol')
disp(' ')
disp('Nilai Eigen Dari Sistem :')
Eigen_R_C = abs(eig(Ar_C))
%inisiasi, jumlah eigen stabil dan tidak stabil
Tidak_Stabil_R_C = 0;
Stabil_R_C = 0;
Stabil_Asimptotik_R_C = 0;
m = size(Ar_C);
for i = 1:m
    if real(Eigen_R_C(i)) > 1
       Tidak_Stabil_R_C = Tidak_Stabil_R_C +1;
    end
        if real(Eigen_R_C(i)) < 1
            Stabil_Asimptotik_R_C = Stabil_Asimptotik_R_C +1;
        end
            if real(Eigen_R_C(i)) == 1
               Stabil_R_C = Stabil_R_C +1;
            end
end
Stabil_R_C = Stabil_R_C + Stabil_Asimptotik_R_C;
Tidak_Stabil_R_C;
Stabil_R_C;

if Stabil_R_C == m
    disp('Maka Sistem Stabil')
else
    disp('Maka Sistem Tidak Stabil')
end
disp(' ')
disp('=========================================')
fprintf('Uji Keterkontrolan dan Uji Keteramatan Sistem Terkontrol \n\n')
%Keterkontrolan
Mc_R_C = ctrb(sys_reduksi_kontrol);
RankMc_R_C = rank(Mc_R_C);
%Ketermatan
Mo_R_C = obsv(sys_reduksi_kontrol);
RankMo_R_C = rank(Mo_R_C);

if rank(Ar_C)==rank(Mc_R_C)
    if rank(Ar_C)==rank(Mo_R_C)
        fprintf('rank A = rank Mc = rank Mo = %.0f\n',n');
        fprintf('Maka Sistem Tereduksi Dengan Kontrol, terkontrol dan teramati\n');
    else
        fprintf('rank A = rank Mc, rank A =/= rank Mo\n\n');
        fprintf('Maka Sistem Tereduksi Dengan Kontrol, terkontrol namun tidak teramati\n');
    end
else rank(A)~= rank(Mc)
    if rank(A)==rank(Mo)
       fprintf('rank A =/= rank Mc, rank A = rank Mo\n\n');
       fprintf('Maka Sistem Tereduksi Dengan Kontrol, tidak terkontrol namun teramati\n');
    else
       fprintf('rank A =/= rank Mc =/= rank Mo\n\n');
       fprintf('Maka Sistem Tereduki Dengan Kontrol, tidak terkontrol dan tidak teramati\n');
    end
end

%'=========================================')
%Gramian Keterkontrolan
W_R_C = gram(sys_reduksi_kontrol,'c');
%Cek W definit positif
EigenW_R_C = eig(W_R_C);
%Cek W nonsingular
DeterminanW_R_C = det(W_R_C);
%'=========================================')
%Gramian Keteramatan
M_R_C = gram(sys_reduksi_kontrol,'o');
%Cek M definit positif
EigenW_R_C = eig(M_R_C);
%Cek M nonsingular
DeterminanW_R_C = det(M_R_C);
disp('')

toc
 %% Perhitungan Error
disp(' ========================================= ')
disp (' ')
disp(' Nilai Error Sistem Tereduksi Yang Dikontrol ')
Error2 = norm(G_awal - G_Tereduksi_Kontrol)

%% Grafik

% Grafik Frekuensi Respon Sistem Dengan Kontrol dan Sistem Dengan Kontrol Yang Tereduksi

figure(1);
t=logspace(-3,3,200);
[mag,pha]=bode(A,B,C,D,1,t); %Sistem Awal
[magc,phac]=bode(Ac,Bc,Cc,Dc,1,t); %Sistem dengan Kontrol
[magrs,phars]=bode(Ar_C,Br_C,Cr_C,Dr_C,1,t);%Sistem tereduksi yang dikontrol

    if Norm_BT <= sigma_r %Reduksi -> Kontrol
       semilogx(t,20*log10(mag),'y-',t,20*log10(magc),'k-.',t,20*log10(magrs),'b-.')
       legend('Sistem Awal','Sistem Awal Dengan Kontrol','Sistem Tereduksi Yang Dikontrol')
    else
        semilogx(t,20*log10(mag),'y-',t,20*log10(magc),'k-.')
        legend('Sistem Awal','Sistem Awal Dengan Kontrol')
    end 
title(['Frekuensi Response Sistem Awal dan Sistem Dengan Orde ' num2str(r)]); 
xlabel('Frequency')
ylabel('Gain')
hold on

figure(2)
sys_1 = ss(A,B,C,D,1);%Sistem Awal
sys_2 = ss(Ac,Bc,Cc,Dc,1);%Sistem dengan Kontrol
sys_R_C = ss(Ar_C,Br_C,Cr_C,Dr_C,1); %Sistem tereduksi yang dikontrol
stepplot(sys_1,'r-',sys_2,'g-',sys_R_C,'y-.')
legend('Sistem Awal','Sistem Dengan Kontrol','Sistem Tereduksi Yang Dikontrol')


figure(3);
w=logspace(0,2,500);
P=pck(A,B,C,D); %Sistem Awal
Pc=pck(Ac,Bc,Cc,Dc); %Sistem Dengan Kontrol
Prc=pck(Ar_C,Br_C,Cr_C,Dr_C); %Sistem tereduksi yang dikontrol
Ger=msub(P,Pc);
Gerrc=msub(P,Prc);
Gf1=frsp(Ger,w);
Gf3=frsp(Gerrc,w);
[u1,s1,v1]=vsvd(Gf1);
[u3,s3,v3]=vsvd(Gf3);
set(0,'defaultfigurecolor',[1 1 1])

    if Norm_BT <= sigma_r %Reduksi -> Kontrol
        vplot('liv,m',s1,s3,':');
        legend('Sistem Dengan Kontrol','Sistem Tereduksi Yang Dikontrol')
    else
        vplot('liv,m',s1,':');
        legend('Sistem Dengan Kontrol')
    end 

title(['Error sistem orde ' num2str(r)])
xlabel('Frequency')
ylabel('Error')