#!/usr/bin/octave -f
arg_list = argv();
file_name = arg_list{1};
cmStart = str2double(arg_list{2});
cmEnd = str2double(arg_list{3});
cbStart = str2double(arg_list{4});
cbEnd = str2double(arg_list{5});
nPoints  = str2num(arg_list{6});
m = csvread(strcat('statistics/', file_name, '.csv'));
x = linspace(cmStart, cmEnd, nPoints);
y = linspace(cbStart, cbEnd, nPoints);
imagesc(x, y, m);
set(gca,'YDir','normal');
xlabel('cm');
ylabel('cb');
title('Number of non Horn clauses');
colormap(gray);
colorbar;
print(strcat('plots/', file_name, '.png'), '-dpng');