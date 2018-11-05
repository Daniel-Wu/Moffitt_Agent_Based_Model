all_p(end, :);

time; %??

disp('Creating figure')

figure(2)
subplot(3, 1, 1)
hold on


if (constantAbi == 1)
    abiStart = find(treatmentIndex == 3, 1, 'first');
    plot([abiStart abiStart], [0 300], 'k', 'LineWidth', 4);
else
    for i = 1:1:time
        if ( treatmentIndex(i) == 3 )
            line([i i], [0 300], 'Color', [0.9, 0.9, 0.9])
        end
    end
end

plot(all_x(:,1), 'r' , 'LineWidth', 4, 'DisplayName', 'T+');
plot(all_x(:,2), ':g' , 'LineWidth', 4, 'DisplayName', 'TP');
plot(all_x(:,3), '--b' , 'LineWidth', 4, 'DisplayName', 'T-');
xlabel('Generations', 'FontSize', 16)
ylabel('Population Size', 'FontSize', 16)
% ylim([0 300])
xlim([0 time])
xlim([0 length(all_PSA)])
%legend ({'T+' 'TP' 'T-'})

subplot(3, 1, 2)
hold on

if (constantAbi == 1)
    abiStart = find(treatmentIndex == 3, 1, 'first');
    plot([abiStart abiStart], [0 300],  'k', 'LineWidth', 4);
else
    for i = 1:1:time
        if ( treatmentIndex(i) == 3 )
            line([i i], [0 300], 'Color', [0.9, 0.9, 0.9])
        end
    end
end

plot(all_p(:,1), 'r' , 'LineWidth', 4, 'DisplayName', 'T+');
plot(all_p(:,2), ':g' , 'LineWidth', 4, 'DisplayName', 'TP');
plot(all_p(:,3), '--b' , 'LineWidth', 4, 'DisplayName', 'T-');
xlabel('Generations', 'FontSize', 16)
ylabel('Frequencies', 'FontSize', 16)
xlim([1 time-1])
xlim([0 length(all_PSA)])

ylim([0 1]);
%legend ({'T+' 'TP' 'T-'})


subplot(3, 1, 3)
hold on

if (constantAbi == 1)
    abiStart = find(treatmentIndex == 3, 1, 'first');
    plot([abiStart abiStart], [0 300],  'k', 'LineWidth', 4);
else
    for i = 1:1:time
        if ( treatmentIndex(i) == 3 )
            line([i i], [0 300], 'Color', [0.9, 0.9, 0.9])
        end
    end
end
plot(1:time, all_PSA, 'k' , 'LineWidth', 4);
xlim([0 time])
xlim([0 length(all_PSA)])

ylim([0 max(all_PSA) + 10])
xlabel('Generations', 'FontSize', 16)
ylabel('PSA', 'FontSize', 16)