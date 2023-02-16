


figure;
hold on
plot([freq ],[cr(:,1)],'r+','MarkerSize',5,'LineWidth',2)

plot([freq freq],[cr_t(:,1);cr_t(:,2)],'k-','MarkerSize',5,'LineWidth',1)

xlim([freq(1) freq(end)]);
box on;
set(gca,'LineWidth',1)
xlabel('Frequency (Hz)','Fontname', ' Arial ')
ylabel('Phase velocity (m/s)','Fontname', ' Arial ')
h_vs=legend('Approximate dispersion curves','Theoretical leaky modes');
