clc;
clear all;

disp('========== Simulasi Sistem Konduksi Panas ==========')
disp('Tahapan 1 : Sistem Diberikan Kontrol Kemudian Direduksi ')
disp(' ')

tic

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
C(1,2)=1;
C(1,3)=1;
C(1,4)=1;
%Kontruksi matriks D
D=[1];
disp('')
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
disp('')
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
disp('')

%% Desain Kontrol Dengan LQR

disp(' ')
disp('Desain Kontrol LQR')
%Kontrol Sistem Awal
Q = C'*C;
R = 1;
[K] = dlqr(A,B,Q,R)


%% Sistem Terkontrol
disp('=========================================')
disp('')
disp('========== Sistem Terkontrol ==========')
Ac = A-B*K;
Bc = B;
Cc = C;
Dc = D;
sys_terkontrol = ss(Ac,Bc,Cc,Dc,1)

% Fungsi transfer sistem awal
G_Sc = tf(sys_terkontrol); 

%% Analsisa Sistem Dengan Desain Kontrol

disp('=========================================')
disp('Analisa Sistem Terkontrol')
disp(' ')
disp('Nilai Eigen Dari Sistem :')
Eigen_Sc = abs(eig(Ac))
%inisiasi, jumlah eigen stabil dan tidak stabil
Tidak_Stabil_Sc = 0;
Stabil_Sc = 0;
Stabil_Asimptotik_Sc = 0;
m = size(Ac);
for i = 1:m
    if real(Eigen_Sc(i)) > 1
       Tidak_Stabil_Sc = Tidak_Stabil_Sc +1;
    end
        if real(Eigen_Sc(i)) < 1
            Stabil_Asimptotik_Sc = Stabil_Asimptotik_Sc +1;
        end
            if real(Eigen_Sc(i)) == 1
               Stabil_Sc = Stabil_Sc +1;
            end
end
Stabil_Sc = Stabil_Sc + Stabil_Asimptotik_Sc;
Tidak_Stabil_Sc;
Stabil_Sc;

if Stabil_Sc == m
    disp('Maka Sistem Stabil')
else
    disp('Maka Sistem Tidak Stabil')
end
disp(' ')
disp('=========================================')
fprintf('Uji Keterkontrolan dan Uji Keteramatan Sistem Terkontrol \n\n')
%Keterkontrolan
Mc_Sc = ctrb(sys_terkontrol);
RankMc_Sc = rank(Mc_Sc);
%Ketermatan
Mo_Sc = obsv(sys_terkontrol);
RankMo_Sc = rank(Mo_Sc);

if rank(A)==rank(Mc)
    if rank(A)==rank(Mo)
        fprintf('rank A = rank Mc = rank Mo = %.0f\n',n');
        fprintf('Maka Sistem Awal Dengan Kontrol, terkontrol dan teramati\n');
    else
        fprintf('rank A = rank Mc, rank A =/= rank Mo\n\n');
        fprintf('Maka Sistem Awal Dengan Kontrol, terkontrol namun tidak teramati\n');
    end
else rank(A)~= rank(Mc)
    if rank(A)==rank(Mo)
       fprintf('rank A =/= rank Mc, rank A = rank Mo\n\n');
       fprintf('Maka Sistem Awal Dengan Kontrol, tidak terkontrol namun teramati\n');
    else
       fprintf('rank A =/= rank Mc =/= rank Mo\n\n');
       fprintf('Maka Sistem Awal Dengan Kontrol, tidak terkontrol dan tidak teramati\n');
    end
end

%'=========================================')
%Gramian Keterkontrolan
W_Sc = gram(sys_terkontrol,'c');
%Cek W definit positif
EigenW_Sc = eig(W_Sc);
%Cek W nonsingular
DeterminanW_Sc = det(W_Sc);
%'=========================================')
%Gramian Keteramatan
M_Sc = gram(sys_terkontrol,'o');
%Cek M definit positif
EigenW_Sc = eig(M_Sc);
%Cek M nonsingular
DeterminanW_Sc = det(M_Sc);
disp('')


%% Kontruksi matriks T 

%Menentukan matriks phi dimana berlaku W=phi(transpose)*phi'
phi_Sc=chol(W_Sc);

%Cek W=phi(transpose)*phi
Cek_Sc = phi_Sc'*phi_Sc;

%Diagonalisai phi*M*phi(transpose) sedemikian hingga berlaku phi*M*phi(transpose) = U*(sigma^2)*U(transpose)
J_Sc = phi_Sc*M_Sc*phi_Sc';
[U_Sc,S_Sc,V_Sc]=svd(J_Sc);

%muncul J_Sc=U*S*V' dan S_Sc = sigma^2
%yang dicari diagonalisai J sehingga J=U*(sigma^2)*U'
%sehingga akar S menunjukkan nilai singular hankel
sigma_Sc = (S_Sc).^(1/2);

%Diperoleh Matriks Transformasi T
T_Sc=phi_Sc'*U_Sc*(sigma_Sc)^(-0.5);

%% Kontruksi Matriks Setimbang

disp(' ')
disp('========== Sistem Terkontrol Yang Setimbang ==========')
disp(' ')
As_Sc=(inv(T_Sc))*Ac*T_Sc;
Bs_Sc=(inv(T_Sc))*Bc;
Cs_Sc=Cc*T_Sc;
Ds_Sc=Dc;
sys_setimbang_Sc=ss(As_Sc,Bs_Sc,Cs_Sc,Ds_Sc,1)

%% Analisa Sistem Terkontrol yang Setimbang

%Fungsi transfer sistem setimbang
G_setimbang_Sc = tf(sys_setimbang_Sc);

disp('=========================================')
disp('Cek Kestabilan Sistem Terkontrol Yang Setimbang')
Eigen_Setimbang_Sc = abs(eig(As_Sc))
%inisiasi, jumlah eigen stabil dan tidak stabil
Tidak_Stabil_Setimbang_Sc = 0;
Stabil_Setimbang_Sc = 0;
Stabil_Asimptotik_Setimbang_Sc = 0;
n = size(As_Sc);
for i = 1:n
    if real(Eigen_Setimbang_Sc(i)) > 1
       Tidak_Stabil_Setimbang_Sc = Tidak_Stabil_Setimbang_Sc +1;
    end
        if real(Eigen_Setimbang_Sc(i)) < 1
            Stabil_Asimptotik_Setimbang_Sc = Stabil_Asimptotik_Setimbang_Sc +1;
        end
            if real(Eigen_Setimbang_Sc(i)) == 1
               Stabil_Setimbang_Sc = Stabil_Setimbang_Sc +1;
            end
end
Stabil_Setimbang_Sc = Stabil_Setimbang_Sc + Stabil_Asimptotik_Setimbang_Sc;
Tidak_Stabil_Setimbang_Sc;
Stabil_Setimbang_Sc;

if Stabil_Setimbang_Sc == n
    disp('Maka Sistem Stabil')
else
    disp('Maka Sistem Tidak Stabil')
end

disp('')
disp('=========================================')
fprintf('Uji Keterkontrolan dan Uji Keteramatan Sistem Setimbang Dengan Kontrol \n\n')
%Keterkontrolan
Mc_s_Sc = ctrb(sys_setimbang_Sc);
RankMc_s_Sc = rank(Mc_s_Sc);
%Ketermatan
Mo_s_Sc = obsv(sys_setimbang_Sc);
RankMo_s_Sc = rank(Mo_s_Sc);

if rank(A)==rank(Mc)
    if rank(A)==rank(Mo)
        fprintf('rank A  = rank Mc = rank Mo = %.0f\n',n');
        fprintf('Maka Sistem Terkontrol Yang Setimbang, terkontrol dan teramati\n');
    else
        fprintf('rank A  = rank Mc  , rank A  =/= rank Mo  \n\n');
        fprintf('Maka Sistem Terkontrol Yang Setimbang, terkontrol namun tidak teramati\n');
    end
else rank(A)~= rank(Mc)
    if rank(A)==rank(Mo)
       fprintf('rank A  =/= rank Mc  , rank A = rank Mo  \n\n');
       fprintf('Maka Sistem Terkontrol Yang Setimbang, tidak terkontrol namun teramati\n');
    else
      fprintf('rank A  =/= rank Mc  =/= rank Mo  \n\n');
      fprintf('Maka Sistem Terkontrol Yang Setimbang, tidak terkontrol dan tidak teramati\n');
    end
end
disp('')
disp('=========================================')
disp('Gramian Keterkontrolan Sistem Setimbang')
Ws_Sc = gram(sys_setimbang_Sc,'c')

disp('=========================================')
disp('Gramian Keteramatan Sistem Setimbang')
Ms_Sc = gram(sys_setimbang_Sc,'o')

disp('=========================================')
disp('Nilai Singular Hankel Sistem Setimbang')
hsv_Sc = hsvd(sys_setimbang_Sc);
disp(hsv_Sc)

%% Reduksi Sistem Setimbang Dengan Kontrol

k = size(As_Sc);
disp(['Masukkan ukuran sistem tereduksi yang diinginkan (ukuran tidak lebih dari orde matriks A)']);
r = input('r = '); 

%partisi matriks sistem setimbang
A11_Sc = As_Sc(1:r,1:r);
A12_Sc = As_Sc(1:r,r+1:k);
A21_Sc = As_Sc(r+1:k,1:r);
A22_Sc = As_Sc(r+1:k,r+1:k);

B1_Sc = Bs_Sc(1:r);
B2_Sc = Bs_Sc(r+1:k);

C1_Sc = Cs_Sc(1:r);
C2_Sc = Cs_Sc(1+r:k);

Sigma1_Sc = sigma_Sc(1:r,1:r);
Sigma2_Sc = sigma_Sc(r+1,r+1);
sigma_r_Sc = 0;
for i=1+r:k
    sigma_k_Sc = hsv_Sc(i);
    sigma_j_Sc = 2*sigma_k_Sc + sigma_r_Sc;
    sigma_r_Sc = sigma_j_Sc;
end

Ac_R = A11_Sc;
Bc_R = B1_Sc;
Cc_R = C1_Sc;
Dc_R = Dc;

disp('========= Sistem Terkontrol Yang Direduksi =========')
Sys_Tereduksi_Sc = ss(Ac_R,Bc_R,Cc_R,Dc_R,1)
G_TereduksiBT_Sc = tf(Sys_Tereduksi_Sc);

%% Analisa Sistem Terkontrol Yang Direduksi
disp('=========================================')
disp('Analisa Sistem Terkontrol Yang Direduksi')

disp(' ')
disp('Nilai Eigen Sistem Terkontrol Yang Direduksi :')
Eigen_tereduksi_Sc = abs(eig(Ac_R))
%inisiasi, jumlah eigen stabil dan tidak stabil
Tidak_Stabil_R_Sc = 0;
Stabil_R_Sc = 0;
Stabil_Asimptotik_R_Sc = 0;
n = size(Ac_R);
for i = 1:n
    if real(Eigen_tereduksi_Sc(i)) > 1
       Tidak_Stabil_R_Sc = Tidak_Stabil_R_Sc +1;
    end
        if real(Eigen_tereduksi_Sc(i)) < 1
            Stabil_Asimptotik_R_Sc = Stabil_Asimptotik_R_Sc +1;
        end
            if real(Eigen_tereduksi_Sc(i)) == 1
               Stabil_R_Sc = Stabil_R_Sc +1;
            end
end
Stabil_R_Sc = Stabil_R_Sc + Stabil_Asimptotik_R_Sc;
Tidak_Stabil_R_Sc;
Stabil_R_Sc;

if Stabil_R_Sc == n
    disp('Maka Sistem Stabil')
else
    disp('Maka Sistem Tidak Stabil')
end

disp('')
disp('=========================================')
fprintf('Uji Keterkontrolan dan Uji Keteramatan Sistem Terkontrol Yang Direduksi \n\n')
%Keterkontrolan
Mc_Tereduksi_Sc = ctrb(Sys_Tereduksi_Sc);
RankMc_Tereduksi_Sc = rank(Mc_Tereduksi_Sc);
%Ketermatan
Mo_Tereduksi_Sc = obsv(Sys_Tereduksi_Sc);
RankMo_Tereduksi_Sc = rank(Mo_Tereduksi_Sc);

if rank(Ac_R)==rank(Mc_Tereduksi_Sc)
    if rank(Ac_R)==rank(Mo_Tereduksi_Sc)
        fprintf('rank A  = rank Mc = rank Mo = %.0f\n',r');
        fprintf('Maka Sistem Terkontrol Yang Direduksi, terkontrol dan teramati\n');
    else
        fprintf('rank A = rank Mc , rank A =/= rank Mo \n\n');
        fprintf('Maka Sistem Terkontrol Yang Direduksi, terkontrol namun tidak teramati\n');
    end
else rank(Ac_R) ~= rank(Mo_Tereduksi_Sc)
    if rank(Ac_R)==rank(Mo_Tereduksi_Sc) 
       fprintf('rank A =/= rank Mc , rank A  = rank Mo \n\n');
       fprintf('Maka Sistem Terkontrol Yang Direduksi, tidak terkontrol namun teramati\n');
    else
       fprintf('rank A Sistem Terkontrol Yang Direduksi =/= rank Mc Sistem Terkontrol Yang Direduksi =/= rank Mo Sistem Terkontrol Yang Direduksi \n\n');
       fprintf('Maka Sistem Terkontrol Yang Direduksi, tidak terkontrol dan tidak teramati\n');
    end
end
disp('')

%Gramian Keterkontrolan
W_Tereduksi_Sc = gram(Sys_Tereduksi_Sc,'c');
%Cek W definit positif
EigenW_Tereduksi_Sc = eig(W_Tereduksi_Sc);
%Cek W nonsingular
DeterminanW_Tereduksi_Sc = det(W_Tereduksi_Sc);

%Gramian Keteramatan
M_Tereduksi_Sc = gram(Sys_Tereduksi_Sc,'o');
%Cek M definit positif
EigenW_Tereduksi_Sc = eig(M_Tereduksi_Sc);
%Cek M nonsingular
DeterminanW_Tereduksi_Sc = det(M_Tereduksi_Sc);

disp('=========================================')

%Cek eror sistem tereduksi
Error_BT_Sc = G_Sc - G_TereduksiBT_Sc;

disp('Cek Syarat Reduksi Dengan Balanced Truncation Norm_BT_Sc =< NilaiSigma_r_Sc')
Norm_BT_Sc=norm(Error_BT_Sc,inf)
NilaiSigma_r_Sc = sigma_r_Sc
if Norm_BT_Sc <= NilaiSigma_r_Sc
    disp('Syarat Terpenuhi')
else
    disp('Syarat Tidak Terpenuhi')
end

disp('')
disp('Nilai error sistem tereduksi dengan Balanced Truncation : ')
disp(Norm_BT_Sc)

toc

%% Tahap 2

disp('=====================================================================================================')
disp(' ')
disp('========== Simulasi Sistem Konduksi Panas ==========')
disp('Tahapan 2 : Sistem Direduksi Kemudian Dikontrol')
disp(' ')

tic

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
C(1,2)=1;
C(1,3)=1;
C(1,4)=1;
%Kontruksi matriks D
D=[1];
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

%% Analisa Sistem Setimbang Konduksi Panas

%Fungsi transfer sistem setimbang
G_setimbang = tf(sys_s);

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
Qs = Cr'*Cr;
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
disp(' ')
disp(' Nilai Error Sistem Terkontrol Yang Direduksi ')
Error1 = norm(G_awal - G_TereduksiBT_Sc)
disp (' ')
disp(' Nilai Error Sistem Tereduksi Yang Dikontrol ')
Error2 = norm(G_awal - G_Tereduksi_Kontrol)
%% Grafik

% Grafik Frekuensi Respon Sistem Dengan Kontrol dan Sistem Dengan Kontrol Yang Tereduksi

figure(1);
t=logspace(-3,3,200);
[mag,pha]=bode(A,B,C,D,1,t); %Sistem Awal
[magc,phac]=bode(Ac,Bc,Cc,Dc,1,t); %Sistem dengan Kontrol
[magrs_Sc,phars_Sc]=bode(Ac_R,Bc_R,Cc_R,Dc_R,1,t);%Sistem terkontrol yang direduksi
[magrs,phars]=bode(Ar_C,Br_C,Cr_C,Dr_C,1,t);%Sistem tereduksi yang dikontrol


if Norm_BT_Sc <= NilaiSigma_r_Sc %Terkontrol -> Reduksi
    if Norm_BT <= sigma_r %Reduksi -> Kontrol
       semilogx(t,20*log10(mag),'y-',t,20*log10(magc),'k-.',t,20*log10(magrs_Sc),'r-',t,20*log10(magrs),'b-.')
       legend('Sistem Awal','Sistem Awal Dengan Kontrol','Sistem Terkontrol Yang Direduksi','Sistem Tereduksi Yang Dikontrol')
    else
        semilogx(t,20*log10(mag),'y-',t,20*log10(magc),'k-.',t,20*log10(magrs_Sc),'r-')
        legend('Sistem Awal','Sistem Awal Dengan Kontrol','Sistem Terkontrol Yang Direduksi')
    end 
else
     if Norm_BT <= sigma_r %Reduksi -> Kontrol
        semilogx(t,20*log10(mag),'y-',t,20*log10(magc),'k-.',t,20*log10(magrs),'b-.')
        legend('Sistem Awal','Sistem Awal Dengan Kontrol','Sistem Tereduksi Yang Dikontrol')
     else
        semilogx(t,20*log10(mag),'y-',t,20*log10(magc),'k-.')
        legend('Sistem Awal','Sistem Awal Dengan Kontrol')
     end     
end
title(['Frekuency Response Sistem Awal dan Sistem Dengan Orde ' num2str(r)]); 
xlabel('Frequency')
ylabel('Gain')
hold on

figure(2)
sys_1 = ss(A,B,C,D,1);%Sistem Awal
sys_2 = ss(Ac,Bc,Cc,Dc,1);%Sistem dengan Kontrol
sys_C_R = ss(Ac_R,Bc_R,Cc_R,Dc_R,1); %Sistem Dengan Kontrol Yang Direduksi
sys_R_C = ss(Ar_C,Br_C,Cr_C,Dr_C,1); %Sistem tereduksi yang dikontrol
stepplot(sys_1,'r-',sys_2,'g-',sys_C_R,'m-.',sys_R_C,'y-.')
legend('Sistem Awal','Sistem Awal Dengan Kontrol','Sistem Terkontrol Yang Direduksi','Sistem Tereduksi Yang Dikontrol')


figure(3);
w=logspace(0,2,500);
P=pck(A,B,C,D); %Sistem Awal
Pc=pck(Ac,Bc,Cc,Dc); %Sistem Dengan Kontrol
Pr=pck(Ac_R,Bc_R,Cc_R,Dc_R); %Sistem terkontrol yang direduksi
Prc=pck(Ar_C,Br_C,Cr_C,Dr_C); %Sistem tereduksi yang dikontrol
Ger=msub(P,Pc);
Gerc=msub(P,Pr);
Gerrc=msub(P,Prc);
Gf1=frsp(Ger,w);
Gf2=frsp(Gerc,w);
Gf3=frsp(Gerrc,w);
[u1,s1,v1]=vsvd(Gf1);
[u2,s2,v2]=vsvd(Gf2);
[u3,s3,v3]=vsvd(Gf3);
set(0,'defaultfigurecolor',[1 1 1])

if Norm_BT_Sc <= NilaiSigma_r_Sc %Terkontrol -> Reduksi
    if Norm_BT <= sigma_r %Reduksi -> Kontrol
        vplot('liv,m',s1,s2,s3,':');
        legend('Sistem Dengan Kontrol','Sistem Terkontrol Yang Direduksi','Sistem Tereduksi Yang Dikontrol')
    else
        vplot('liv,m',s1,s2,':');
        legend('Sistem Dengan Kontrol','Sistem Terkontrol Yang Direduksi')
    end 
else
     if Norm_BT <= sigma_r %Reduksi -> Kontrol
        vplot('liv,m',s1,s3,':');
        legend('Sistem Dengan Kontrol','Sistem Tereduksi Yang Dikontrol')
     else
        vplot('liv,m',s1,':');
        legend('Sistem Dengan Kontrol')
     end     
end
title(['Error sistem terkontrol yang direduksi menjadi orde ' num2str(r)])
xlabel('Frequency')
ylabel('Error')





