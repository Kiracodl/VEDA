rappture_reader_single('foo.xml')
%%


figureset(4, 3);

plot( Hsample_x, Hsample_y,'r', Topog_x, Topog_y,'b', 'LineWidth', 2);
xlabel('X distance (nm)');
ylabel('Height (nm)');
legend('Actual sample', 'Reported Topography')
grid on;

print -dpdf AMS_Ex3_topo2


figureset(4, 3);

plot(  MeanForce_x, MeanForce_y * 1000,'b', 'LineWidth', 2);
xlabel('X distance (nm)');
ylabel('Mean Force (pN)');
grid on;
ylim([-30 30])
print -dpdf AMS_Ex3_topo2

