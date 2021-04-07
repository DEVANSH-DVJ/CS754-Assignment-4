clc;
clear;
close all;

rng(0);

dctmat = dctmtx(256);

ss = zeros(100,1);
error1 = zeros(100,1);
error2 = zeros(100,1);

for i=1:100
    s = 25;
    ind1 = randi(256, s, 1);
    coeff1 = zeros(256,1);
    coeff1(ind1) = rand(s,1)*100;
    ind2 = randi(256, s, 1);
    coeff2 = zeros(256,1);
    coeff2(ind2) = rand(s,1)*100;

    f1 = dctmat*coeff1;
    f2 = coeff2;
    f = f1 + f2;

    sigma = 0.01 * i * abs(mean(f));

    A = [dctmat eye(256)];
    x = omp_error(A, f, sigma);
    coeff1_recon = x(1:256);
    coeff2_recon = x(257:512);

    f1_recon = dctmat*coeff1_recon;
    f2_recon = coeff2_recon;

    ss(i) = 0.01 * i;
    error1(i) = norm(f1_recon - f1)/norm(f1);
    error2(i) = norm(f2_recon - f2)/norm(f2);
end

% Correlation plot
figure;
plot(ss, error1, 'b');
xlabel("\sigma");
ylabel("$$||\hat{f_1} - f1||/||f1||$$", 'interpreter', 'latex');
title("Varying \sigma keeping s=25");
saveas(gcf, "plots/error1_sigma.jpg");

% Correlation plot
figure;
plot(ss, error2, 'r');
xlabel("\sigma");
ylabel("$$||\hat{f_2} - f_2||/||f_2||$$", 'interpreter', 'latex');
title("Varying \sigma keeping s=25");
saveas(gcf, "plots/error2_sigma.jpg");














ss = zeros(100,1);
error1 = zeros(100,1);
error2 = zeros(100,1);

for i=1:100
    s = 1*i;
    ind1 = randi(256, s, 1);
    coeff1 = zeros(256,1);
    coeff1(ind1) = rand(s,1)*100;
    ind2 = randi(256, s, 1);
    coeff2 = zeros(256,1);
    coeff2(ind2) = rand(s,1)*100;

    f1 = dctmat*coeff1;
    f2 = coeff2;
    f = f1 + f2;

    sigma = 0.01 * abs(mean(f));

    A = [dctmat eye(256)];
    x = omp_error(A, f, sigma);
    coeff1_recon = x(1:256);
    coeff2_recon = x(257:512);

    f1_recon = dctmat*coeff1_recon;
    f2_recon = coeff2_recon;

    ss(i) = s;
    error1(i) = norm(f1_recon - f1)/norm(f1);
    error2(i) = norm(f2_recon - f2)/norm(f2);
end

% Correlation plot
figure
plot(ss, error1, 'b');
xlabel("sparsity level");
ylabel("$$||\hat{f_1} - f1||/||f1||$$", 'interpreter', 'latex');
title("Varying sparsity level keeping \sigma=0.01*mean(f_1+f_2)");
saveas(gcf, "plots/error1_sparsity.jpg");

% Correlation plot
figure
plot(ss, error2, 'r');
xlabel("sparsity level");
ylabel("$$||\hat{f_2} - f_2||/||f_2||$$", 'interpreter', 'latex');
title("Varying sparsity level keeping \sigma=0.01*mean(f_1+f_2)");
saveas(gcf, "plots/error2_sparsity.jpg");








ss = zeros(100,1);
error1 = zeros(100,1);
error2 = zeros(100,1);

for i=1:100
    s = 25;
    ind1 = randi(256, s, 1);
    coeff1 = zeros(256,1);
    coeff1(ind1) = rand(s,1)*100;
    ind2 = randi(256, s, 1);
    coeff2 = zeros(256,1);
    coeff2(ind2) = rand(s,1)*100*i;

    f1 = dctmat*coeff1;
    f2 = coeff2;
    f = f1 + f2;

    sigma = 0.01 * abs(mean(f));

    A = [dctmat eye(256)];
    x = omp_error(A, f, sigma);
    coeff1_recon = x(1:256);
    coeff2_recon = x(257:512);

    f1_recon = dctmat*coeff1_recon;
    f2_recon = coeff2_recon;

    ss(i) = i;
    error1(i) = norm(f1_recon - f1)/norm(f1);
    error2(i) = norm(f2_recon - f2)/norm(f2);
end

% Correlation plot
figure
plot(ss, error1, 'b');
xlabel("k");
ylabel("$$||\hat{f_1} - f1||/||f1||$$", 'interpreter', 'latex');
title("Varying magnitude of f_2 w.r.t f_1 (ratio = k)");
saveas(gcf, "plots/error1_k.jpg");

% Correlation plot
figure
plot(ss, error2, 'r');
xlabel("k");
ylabel("$$||\hat{f_2} - f_2||/||f_2||$$", 'interpreter', 'latex');
title("Varying magnitude of f_2 w.r.t f_1 (ratio = k)");
saveas(gcf, "plots/error2_k.jpg");
