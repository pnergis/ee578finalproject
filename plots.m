%% Plotting for final project

figure(1)
plot(FOM_TE_1(1:56)*100, 'r', 'Linewidth',2)
title('TE mode, Efficiency vs Iteration, 29.5 min')
grid on
xlabel('Iteration');
ylabel('Field Efficiency (%)');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18)
%%
figure(2)
plot(FOM_TM_2(1:96)*100, 'b', 'Linewidth',2)
title('TM mode, Efficiency vs Iteration, 35.2 min')
grid on
xlabel('Iteration');
ylabel('Field Efficiency (%)');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18)
%%
figure(6)
plot(FOM_TE_3(1:296)*100, 'r', 'Linewidth',2)
hold on
plot(FOM_TM_3(1:296)*100, 'b', 'Linewidth',2)
title('TM and TE mode, Efficiency vs Iteration, 4.5 hrs')
grid on
xlabel('Iteration');
ylabel('Field Efficiency (%)');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18)
legend('TE','TM', 'Location', 'southeast')
%%
figure (3)
pcolor(pbs_1)
xlabel('x')
ylabel('y')
title('\epsilon TE');
set(gca,'YDir','normal')
%%
figure (4)
pcolor(pbs_2)
xlabel('x')
ylabel('y')
title('\epsilon TM');
set(gca,'YDir','normal')
%%
figure (5)
pcolor(pbs_3)
xlabel('x')
ylabel('y')
title('\epsilon TM and TE');
set(gca,'YDir','normal')