T = readmatrix("last_folding_energy.csv");
T_random = readmatrix("folding_energy_random.csv");
SLIDING_WINDOW_SIZE = 40;

max_length = 80;
ORF_SIZE = 550;

mean_vals = mean(T);
std_vals = std(T);

mean_vals_rand = mean(T_random);
std_vals_rand = std(T_random);

figure;
x = linspace(-80,length(mean_vals)-80,length(mean_vals)-40);
plot(x, mean_vals(1:end-40), "LineWidth", 1);
xlim([-80, length(mean_vals)-40]);
hold on;
plot(x, mean_vals_rand(1:end-40), "LineWidth", 1, "Color", "red");
xline(0, 'k--', 'LineWidth', 1)
figure;
plot(x, std_vals(1:end-40));
xlim([-80, length(mean_vals)-40]);
hold on;
plot(x, std_vals_rand(1:end-40), "LineWidth", 1, "Color", "red");
xline(0, 'k--', 'LineWidth', 1)