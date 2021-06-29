% GVR 2.20

% Part a

fprintf('Part a.\n')

F1 = 100;
mx5 = 0.926 * 0.15 * F1;
F5 = mx5 / 0.776;
mt4 = 0.92 * 0.55 * F1;
F4 = mt4 / 0.946;
F2 = F1 - F4 - F5;

xt5 = (55 - 0.0454*F2 - mt4) / F5;
xx4 = (15 - 0.0106*F2 - 0.776*F5) / F4;

xb4 = 1 - mt4/F4 - xx4;
xb5 = 1 - xt5 - mx5/F5;

fprintf('Stream 4: wb4 = %.3f, wt4 = %.3f, wx4 = %.3f\n', xb4, mt4/F4, xx4)
fprintf('Stream 5: wb4 = %.3f, wt4 = %.3f, wx4 = %.3f\n', xb5, xt5, mx5/F5)

% Part b

fprintf('\nPart b.\n')
ratio = (0.944*F2)/(0.3*F1);
fprintf('The percentage recovery is %.3f.\n', ratio)