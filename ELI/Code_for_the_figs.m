% fig 3 is Y fig 4 is X
%Figure 3
ax = gca;
ax.FontSize = 50;
ax.FontName = "Times New Roman";
ax.XLabel.String = "Pixels";
ax.XLabel.FontSize = 60;
ax.YLabel.String = "Pixel response (a.u)";
ax.YLabel.FontSize = 60;
ax.Children(1).LineWidth = 3;
ax.Children(1).LineStyle = "--";
ax.Children(1).DisplayName = "Jinc(y)";
ax.Children(2).LineWidth = 2;
ax.Children(2).DisplayName = "Data Y";
legend


Max_data = max(ax.Children(1).YData);
yv = linspace (13.36-2,19.34-2,101);
y_high = Max_data - 2.89e4;
vert_vec_y = linspace(0,y_high,101);
L_top = line(yv,ones(size(yv))*y_high,"LineWidth" , 3 , "LineStyle", ":");
L_top.HandleVisibility = "callback"; L_top.Color = "red";
L_up = line(yv(1)*ones(size(vert_vec_y)),vert_vec_y,"LineWidth" , 3 , "LineStyle", ":");
L_up.HandleVisibility = "callback"; L_up .Color = "red";

L_down = line(yv(end)*ones(size(vert_vec_y)),vert_vec_y,"LineWidth" , 3 , "LineStyle", ":");
L_down.HandleVisibility = "callback"; L_down .Color = "red";

%Figure 4
ax = gca;
ax.FontSize = 50;
ax.FontName = "Times New Roman";
ax.XLabel.String = "Pixels";
ax.XLabel.FontSize = 60;
ax.YLabel.String = "Pixel response (a.u)";
ax.YLabel.FontSize = 60;
ax.Children(1).LineWidth = 3;
ax.Children(1).LineStyle = "--";
ax.Children(1).DisplayName = "Jinc(x)";
ax.Children(2).LineWidth = 2;
ax.Children(2).DisplayName = "Data X";
legend


Max_data = max(ax.Children(1).YData);
xv =  linspace (11.76-2,19.47-2,101);
x_high = Max_data - 3.05e4;


vert_vec_x = linspace(0,x_high,101);
L_top = line(xv,ones(size(xv))*x_high,"LineWidth" , 3 , "LineStyle", ":");
L_top.HandleVisibility = "callback"; L_top.Color = "red";
L_up = line(xv(1)*ones(size(vert_vec_x)),vert_vec_x,"LineWidth" , 3 , "LineStyle", ":");
L_up.HandleVisibility = "callback"; L_up .Color = "red";

L_down = line(xv(end)*ones(size(vert_vec_x)),vert_vec_x,"LineWidth" , 3 , "LineStyle", ":");
L_down.HandleVisibility = "callback"; L_down .Color = "red";
