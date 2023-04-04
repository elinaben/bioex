T = readmatrix("folding_energy_bioinformatics.csv");
T_random = readmatrix("folding_energy_random_bioinformatics.csv");
SLIDING_WINDOW_SIZE = 40;

UTR_SIZE = -23;
ORF_SIZE = 550;
T(T == 0) = NaN;
T_random(T_random == 0) = NaN;
mean_vals = mean(T, 'omitnan');
std_vals = std(T, 'omitnan');

mean_vals_rand = mean(T_random, 'omitnan');
std_vals_rand = std(T_random, 'omitnan');


% Perform the Wilcoxon rank sum test and get the p-value
[pvalue_utr5, h] = ranksum(mean_vals(1:23), mean_vals_rand(1:23));
[pvalue_orf, h] = ranksum(mean_vals(24:end), mean_vals_rand(24:end));

figure;
x = linspace(UTR_SIZE,length(mean_vals),length(mean_vals));
plot(x, mean_vals(1:end), "LineWidth", 1);
xlim([UTR_SIZE, length(mean_vals)-40]);
hold on;
plot(x, mean_vals_rand(1:end), "LineWidth", 1, "Color", "red");
xline(0, 'k--', 'LineWidth', 1)


left_label = sprintf('p=%.4f', pvalue_utr5);
right_label = sprintf('p=%d', pvalue_orf);
text(-3, -5, left_label, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Color','red','FontSize',10)
text(12, -5.5, right_label, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','FontSize',10)

figure;
plot(x, std_vals(1:end));
xlim([UTR_SIZE, length(mean_vals)-40]);
hold on;
plot(x, std_vals_rand(1:end), "LineWidth", 1, "Color", "red");
xline(0, 'k--', 'LineWidth', 1)

