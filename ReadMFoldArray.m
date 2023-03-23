T = readmatrix("folding_energy.csv");
SLIDING_WINDOW_SIZE = 40;
max_length = SLIDING_WINDOW_SIZE;
ORF_SIZE = 500;

mean_vals = mean(T);
std_vals = std(T);

figure;
x = linspace(-40,500,500);
plot(x, mean_vals(1:end-40), "LineWidth", 1);
xlim([-40, length(mean_vals)-40]);
hold on;
xline(0, 'k--', 'LineWidth', 1)
figure;
plot(x, std_vals(1:end-40));
xlim([-40, length(mean_vals)-40]);
hold on;
xline(0, 'k--', 'LineWidth', 1)