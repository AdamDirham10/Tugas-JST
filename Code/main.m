disp('Algoritma Levenberg-Marquardt / Levenberg-Marquardt-Fletcher');
disp('Contoh: Memperkirakan penjualan pada periode berikutnya berdasarkan data penjualan pada periode sebelumnya');
disp('Diasumsikan ada 4 data barang yang sudah diketahui hasil penjualannya');
disp('Masing-masing barang memiliki data penjualan pada 5 periode sebelumnya');
disp('Maka tentukan prediksi penjualan untuk 3 periode berikutnya');
disp('Diasumsikan data penjualan adalah sebagai berikut: ');
disp('Nama Barang, Periode1, Periode2, Periode3, Periode4, Periode5');
disp('Pensil     , 10      , 12      , 8       , 11      , 20');
disp('Penghapus  , 15      , 12      , 13      , 14      , 14');
disp('Pena       , 9       , 7       , 15      , 15      , 10');
disp('Penggaris  , 1       , 2       , 4       , 12      , 13');
disp(char(10));

namaBarang = {'Pensil' 'Penghapus' 'Pena' 'Penggaris'};
data=[10 12 8 11 20; 15 12 13 14 14; 9 7 15 15 10; 1 2 4 12 13];
jumlahData = size(data,2);

fx = @(x) sqrt(-1*x(1)*x(1) + -0.5*x(2)*x(2) + 0.5*x(3)*x(3) + 0.75*x(4)*x(4) + x(5)*x(5));

epsX  = 1e-7;
epsFX = 1e-7;

jacobian = 'finiteJacobian';

maksIterasi = 500;

minBobotSelisih  = 0.25;
maksBobotSelisih = 0.75;

for i=1:size(data,1)

    output = LMF(fx,data(i,:)', epsX, epsFX, jacobian, maksIterasi, minBobotSelisih, maksBobotSelisih);
    
    lmf(i) = round(output(end),0);
end

disp('Hasil Perhitungan dengan metode Levenberg-Marquardt / Levenberg-Marquardt-Fletcher')
disp('Nama Barang, Periode1, Periode2, Periode3, Periode4, Periode5, Prediksi Periode berikutnya')
for i = 1:size(lmf, 2),
    disp([char(namaBarang(i)), blanks(11 - cellfun('length',namaBarang(i))), ', ', ... 
         num2str(data(i,1)), blanks(8 - length(num2str(data(i,1)))), ', ', ...
         num2str(data(i,2)), blanks(8 - length(num2str(data(i,2)))), ', ', ...
         num2str(data(i,3)), blanks(8 - length(num2str(data(i,3)))), ', ', ...
         num2str(data(i,4)), blanks(8 - length(num2str(data(i,4)))), ', ', ...
         num2str(data(i,5)), blanks(8 - length(num2str(data(i,5)))), ', ', ...
         num2str(lmf(i))])
end