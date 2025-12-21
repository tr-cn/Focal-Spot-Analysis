close all; clear; clc;
global Pixel_size effective_magnification
Plot = true;
Laser_eng = 560;%  [J]
Laser_wavelenth = 0.840; %[um]
pulse_duration = 35*1e-15; %[sec]
window_size = 80;
Directory=  'P:\ELI\Alignment\New folder';
% Directory= 'C:\Users\TamirCo\Desktop\before_experiment\before_experiment';
% Directory = dir('C:\Tamir Cohen Files\תואר שלישי\פרוייקט פמטו\Focuse analizing\90 degrees off axis Focuse Test');

Dir = dir();
f_x = {}; f_y = {};
for ii =3:length(Dir)
    i=ii-2;

    FileName = [Dir(ii).name];

    if ~strcmp(FileName(end-3:end),'tiff')
        continue;
    end
    % FileName =["C:\Users\TamirCo\Desktop\relativeDifferencesBetweenFoci.png"];
    [Image] = figure_loader(FileName,i);
    fig = gcf;
    fig_num = fig.Number;
    % if i == 1
    % [Background] = background_Remover (Image);
    Background = mean(mean(Image([1:150,430:end],:)));
    % end

    Image  = Image - Background;
    Image_full = Image;
    % [Image] = Region_of_intrest (Image);
    % SmoothImag = smoothdata(Image,2,"movmedian");
    [xx,yy] = find(Image == max(max(Image)));
    Image = Image(xx-window_size:xx+window_size,yy-window_size:yy+window_size );
    imagesc(Image); axis equal;
    % Image = Iamge_Smoother (Image);

    old_length = length(Image);
    % Image = interp2(Image,[2,2]);
    % Image = smoothdata(Image,2,"gaussian",[5,5]);

    New_length = length(Image);
    imagesc(Image); axis equal;


    % Pixel_size = 3.45 * old_length/New_length; % um
    magnification = 15.042;
    Pixs_per_mu = 4.36;
    binning = 1;
    Pixel_size = 1/Pixs_per_mu*magnification*binning;
    effective_magnification = magnification;%* 0.9287;

    [Max_y_index,Max_x_index] = find(Image == max(max(Image)));
    Max_x_index = round(mean(Max_x_index));
    Max_y_index = round(mean(Max_y_index));
    % Mean_diff = median(median(Image));
    % vec_x = 1:1920;
    % vec_y = 1:1080;
    % effetive_field = 150;
    %
    % if Max_y_index-effetive_field>=0
    %     vec_y = Max_y_index-effetive_field:Max_y_index+effetive_field;
    % else
    %     vec_y = 0:Max_y_index+effetive_field+abs(effetive_field-Max_y_index);
    % end
    %
    % if Max_x_index-effetive_field>=0
    %     vec_x = Max_x_index-effetive_field:Max_x_index+effetive_field;
    % else
    %     vec_x = 0:Max_x_index+effetive_field+abs(effetive_field-Max_x_index);
    % end
    %
    % Image = Image(vec_y,vec_x);
    %
    % imagesc()
    %
    %
    %
    %
    % I_x = Image(effetive_field+1,:);
    % I_y = Image(:,effetive_field+1)';
    I_x = Image(:,Max_x_index);
    I_y = Image(Max_y_index,:);

    vec_x = 1:size(Image,1) ;
    vec_y = 1:size(Image,2);

    x0 = [max(I_x),Max_x_index,1,0];
    y0 =  [max(I_y),Max_y_index,1,0];



    [center_x(i), f_x{i},f_cont_x{i},pix_vec_fit_x] = bessel_fit (vec_x, I_x,x0);
    [center_y(i), f_y{i},f_cont_y{i},pix_vec_fit_y] = bessel_fit (vec_y, I_y,y0);



    x_microns{i} = To_microns(vec_x);
    y_microns{i} = To_microns(vec_y);

    center_x_microns{i}(i) = To_microns(center_x(i));
    center_y_microns{i}(i) = To_microns(center_y(i));

    % [fwhm_x(i)] = FWHM_fun(x_microns{i},f_x);
    % [fwhm_y(i)] = FWHM_fun(y_microns{i},f_y);
    [fwhm_x_mu(i),fwhm_x_pxl(i),First_zero_x_mu(i),First_zero_x_pxl(i),e_minus2_x_mu(i),e_minus2_x_pxl(i)] = FWHM_and_First_zero_fun_1D(x_microns{i},f_x{i});
    [fwhm_y_mu(i),fwhm_y_pxl(i),First_zero_y_mu(i),First_zero_y_pxl(i),e_minus2_y_mu(i),e_minus2_y_pxl(i)] = FWHM_and_First_zero_fun_1D(y_microns{i},f_y{i});

    C = sqrt(center_x(i)^2 + center_y(i)^2);
    R_FWHM_mu(i)  = sqrt( (0.5*fwhm_x_mu(i))*(0.5*fwhm_y_mu(i)));
    R_FirstZero_mu(i) =  sqrt( (First_zero_x_mu(i)*0.5)*(0.5*First_zero_y_mu(i)));
    R_e_minus2_mu(i)  =  sqrt( (e_minus2_x_mu(i)*0.5)*(0.5*e_minus2_y_mu(i)));
    theta = linspace(0,2*pi,1001);

    Eng_Pers_FWHM(i) = sum_pixels_in_ellipse(Image, [center_x(i),center_y(i)], [fwhm_x_pxl(i)/2,fwhm_y_pxl(i)/2])/sum(sum(Image_full))*100;
    Eng_Pers_e_minus_2(i) = sum_pixels_in_ellipse(Image, [center_x(i),center_y(i)], [e_minus2_x_pxl(i)/2,e_minus2_y_pxl(i)/2])/sum(sum(Image_full))*100;
    Eng_Pers_First_zero(i) = sum_pixels_in_ellipse(Image, [center_x(i),center_y(i)], [First_zero_x_pxl(i)/2,First_zero_y_pxl(i)/2])/sum(sum(Image_full))*100;

    Rad_70_persent_pxl = Seventy_pers_radius (Image,Image_full,[center_x(i),center_y(i)]);
    Dia_70_persent_um = To_microns(Rad_70_persent_pxl) * 2;
    if Plot
        figure(fig_num);
        imagesc(y_microns{i},x_microns{i},Image); axis equal
        title ( Dir(ii).name)
        line(center_y_microns{i}(i)*ones(1,length(x_microns{i})),x_microns{i},'LineWidth',6,'color','white');
        line(y_microns{i},center_x_microns{i}(i)*ones(1,length(y_microns{i})),'LineWidth',6,'color','white');



        % line(center_y_microns{i}(i) + R_FWHM*cos(theta),center_x_microns{i}(i) + R_FWHM*sin(theta),'LineStyle','--','LineWidth',5,'Color','magenta')
        % line(center_y_microns{i}(i) + R_FirstZero*cos(theta),center_x_microns{i}(i) + R_FirstZero*sin(theta),'LineStyle','--','LineWidth',5,'Color','Yellow')
        line(center_y_microns{i}(i) + (0.5*fwhm_y_mu(i))*cos(theta),center_x_microns{i}(i) + (0.5*fwhm_x_mu(i))*sin(theta),'LineStyle','--','LineWidth',5,'Color','magenta')
        line(center_y_microns{i}(i) + (0.5*e_minus2_y_mu(i))*cos(theta),center_x_microns{i}(i) + (0.5*e_minus2_x_mu(i))*sin(theta),'LineStyle','--','LineWidth',5,'Color','green')
        % line(center_y_microns{i}(i) + (0.5*First_zero_y_mu(i))*cos(theta),center_x_microns{i}(i) + (0.5*First_zero_x_mu(i))*sin(theta),'LineStyle','--','LineWidth',5,'Color','Yellow')
        line(center_y_microns{i}(i) + (0.5*Dia_70_persent_um)*cos(theta),center_x_microns{i}(i) + (0.5*Dia_70_persent_um)*sin(theta),'LineStyle','--','LineWidth',5,'Color','Red')

        figure(3)
        plot(vec_x, I_x, pix_vec_fit_x, f_cont_x{i})
        figure(4)
        plot(vec_y, I_y, pix_vec_fit_y, f_cont_y{i})
        figure(5)
        plot(x_microns{i}, f_x{i})
        figure(6)
        plot(y_microns{i}, f_y{i})

    end


    I_fwhm = Laser_eng*(Eng_Pers_FWHM/100)./(pulse_duration.* (fwhm_x_mu./2).*(fwhm_y_mu/2)*(1e-4)^2*pi);
    a0_fwhm = 0.855*1e-9*Laser_wavelenth*sqrt(I_fwhm);

    I_e_minus2 = Laser_eng*(Eng_Pers_e_minus_2/100)./(pulse_duration.* (e_minus2_x_mu./2).*(e_minus2_y_mu/2)*(1e-4)^2*pi);
    a0_e_minus2 = 0.855*1e-9*Laser_wavelenth*sqrt(I_e_minus2);

    I_First_zero = Laser_eng*(Eng_Pers_First_zero/100)./(pulse_duration.* (First_zero_x_mu./2).*(First_zero_y_mu/2)*(1e-4)^2*pi);
    a0_First_zero = 0.855*1e-9*Laser_wavelenth*sqrt(I_First_zero);


    I_70_persent = Laser_eng*(70/100)./(pulse_duration.* (Dia_70_persent_um./2*1e-4)^2*pi);
    a0_70_persent = 0.855*1e-9*Laser_wavelenth*sqrt(I_70_persent);
    
    w0_70_um = Dia_70_persent_um/2 / 0.7759; % 1-exp(2*r^2/w0^2) = 0.7
    I_w0_70 = Laser_eng./((w0_70_um*1e-4).^2*pi*pulse_duration);
    a0_w0_70 = 0.855*1e-9*Laser_wavelenth*sqrt(I_w0_70);








    disp(["a0 - FWHM", a0_fwhm])
    disp (["a0 - e_minus2", a0_e_minus2])
    disp (["a0 - First_zero", a0_First_zero])
    disp (["a0 - 70% of total energy", a0_70_persent])
    disp (["a0 - w0_70% of total energy", a0_w0_70])
end

mean_x_center = mean(center_x);
std_x_center = To_microns(std(center_x));
mean_y_center = mean(center_y);
std_y_center = To_microns(std(center_y));

mean_x_fwhm = mean(fwhm_x_mu);
std_x_fwhm = std(fwhm_x_mu);
mean_y_fwhm = mean(fwhm_y_mu);
std_y_fwhm = std(fwhm_y_mu);

[theta,R_centers] = Radius_of_convergens(mean_x_center, center_x, mean_y_center, center_y);
R_centers = To_microns(R_centers);
figure
polarscatter(theta,R_centers);



function [Image] = figure_loader(FileName,i)
Image = imread(FileName);
% Image = Image(:,:,1) + Image(:,:,2) + Image(:,:,3);
S = size (Image);
if length(S)>2
    if S(3) == 4
        Image = Image(:,:,1:3);
    end
    Image = rgb2gray(Image);
end
Image = double(Image);

% if i == 1
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(Image); axis equal;

% end
end

function rgbImage = bayer2rgb_RGGB(bayerImage)
% bayer2rgb_RGGB Transforms a Bayer pattern image (RGGB) into an RGB image.
%   rgbImage = bayer2rgb_RGGB(bayerImage) takes a grayscale Bayer image
%   and debayers it into a color RGB image using bilinear interpolation.
%   The function assumes an 'RGGB' pattern starting from the top-left pixel.

[rows, cols] = size(bayerImage);

% Initialize the RGB image with zeros
rgbImage = zeros(rows, cols, 3, class(bayerImage));

% Extract individual color components based on the RGGB pattern
% Red pixels
R = zeros(rows, cols, class(bayerImage));
R(1:2:end, 1:2:end) = bayerImage(1:2:end, 1:2:end);

% Green pixels
G = zeros(rows, cols, class(bayerImage));
G(1:2:end, 2:2:end) = bayerImage(1:2:end, 2:2:end);
G(2:2:end, 1:2:end) = bayerImage(2:2:end, 1:2:end);

% Blue pixels
B = zeros(rows, cols, class(bayerImage));
B(2:2:end, 2:2:end) = bayerImage(2:2:end, 2:2:end);

% Interpolate missing values using bilinear interpolation
% Red channel interpolation
R_interp = imfilter(R, [1 2 1; 2 4 2; 1 2 1]/4, 'conv');

% Green channel interpolation
% Green pixels have neighbors of Red and Blue.
% Interpolate Green at Red and Blue locations.
G_interp_at_RB = imfilter(G, [0 1 0; 1 4 1; 0 1 0]/4, 'conv');
% Combine original Green pixels with interpolated ones
G_final = G + G_interp_at_RB; % original G is 0 where interpolation is needed

% Blue channel interpolation
B_interp = imfilter(B, [1 2 1; 2 4 2; 1 2 1]/4, 'conv');

% Combine the interpolated channels into the final RGB image
rgbImage(:,:,1) = R_interp;
rgbImage(:,:,2) = G_final;
rgbImage(:,:,3) = B_interp;

% Handle borders - simple cropping to avoid boundary issues from filtering
rgbImage = rgbImage(2:end-1, 2:end-1, :);

end

function [mean_background] = background_Remover (Image)
msg = msgbox('Mark background area');
uiwait(msg)
rec = drawrectangle("Color",[1 1 1]);
pos_Rect = round(rec.Position);
top_left = pos_Rect(1:2); bottom_right = pos_Rect(1:2) + pos_Rect(3:4); % [X,Y]
X = top_left(1):bottom_right(1);
Y = top_left(2):bottom_right(2);
mean_background = mean(mean(Image(Y,X)));
% Image = Image - mean_background;
end

function [Image_reduced] = Region_of_intrest (Image)
imagesc(Image); axis equal;
msg = msgbox('Mark the region of intrest');
uiwait(msg)
rec = drawrectangle("Color",[1 1 1]);
pos_Rect = round(rec.Position);
top_left = pos_Rect(1:2); bottom_right = pos_Rect(1:2) + pos_Rect(3:4); % [X,Y]
X = top_left(1):bottom_right(1);
Y = top_left(2):bottom_right(2);
Image_reduced = Image(Y,X);


end

function [Image] = Iamge_Smoother (Image)

% Image = smoothdata(Image,'movmedian');
Mean_diff = mean(mean(abs(diff(Image))));
Size = size(Image);
for y = 1:Size(1)
    image_x = Image(y,:);
    diff_x = diff(image_x);
    flag = 0;
    for d = 1:length(diff_x)
        if flag == 1
            flag = 0;
        end

        if ((abs(diff_x(d)) > 100 * Mean_diff) && (flag == 0))
            image_x(d+1) = mean([image_x(d),image_x(d+1)]);
            flag = 1;
        end

    end
    Image(y,:) = image_x;
end

end

function [center,f,f_con,pix_vec] = bessel_fit (pix_vec, I,x0)

pix_vec = reshape(pix_vec,size(I));
Vf = @(x) (pix_vec-x(2)) / x(3);
Bessel_fun = @(x) x(1) *(besselj(1,Vf(x)) ./ Vf(x)).^2 + x(4);
err = @(x) sum((Bessel_fun(x) - I).^2);

options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals',100000);


x = fminsearch(err,x0,options);
% x = fminsearch(err,x0);
center = x(2);
vf = (pix_vec-x(2)) / x(3);
f = x(1) *(besselj(1,vf) ./ vf).^2 + x(4);

pix_vec = linspace(pix_vec(1),pix_vec(end),1001);
vf = (pix_vec-x(2)) / x(3);
f_con = x(1) *(besselj(1,vf) ./ vf).^2 + x(4);

end

function [fwhm] = FWHM_fun(microns,f)

f = f/max(f);
f_MI= find(f==max(f)); % index of the maximum of f
left_side = interp1(f(1:f_MI),microns(1:f_MI),0.5);
right_side = interp1(f(f_MI:end),microns(f_MI:end),0.5);
fwhm = abs(right_side-left_side);
end

function [fwhm_mu,fwhm_pxl,zero_x_mu,zero_x_pxl,e_minus2_mu,e_minus2_pxl] = FWHM_and_First_zero_fun_1D(microns,f)
% f = f-min(f);
% f = f/max(f);
len_f_old = length(f);
min_i_pxl = find(f==min(f));  min_i_pxl=min_i_pxl(1);
max_i_pxl = find(f==max(f));  max_i_pxl=max_i_pxl(1);

% interp_val = 3;
% f = interp(f,interp_val);
len_f_new = length(f);
pxl = 1:length(f);
% microns  = interp(microns,interp_val);
f = f-min(f);
half = max(f)/2;
f_e_minus2 = max(f) * exp(-2);
f_MI= find(f==max(f)); % index of the maximum of f

left_side_fwhm_mu = interp1(f(1:f_MI),microns(1:f_MI),half);
left_side_fwhm_pxl = interp1(f(1:f_MI),1:f_MI,half) * len_f_old/len_f_new;

right_side_fwhm_mu  = interp1(f(f_MI:end),microns(f_MI:end),half);
right_side_fwhm_pxl = interp1(f(f_MI:end),f_MI:length(f),half)* len_f_old/len_f_new;

fwhm_mu  = abs(right_side_fwhm_mu-left_side_fwhm_mu);
fwhm_pxl = abs(right_side_fwhm_pxl-left_side_fwhm_pxl);

left_side_e_minus2_mu  = interp1(f(1:f_MI),microns(1:f_MI),f_e_minus2);
left_side_e_minus2_pxl = interp1(f(1:f_MI),1:f_MI,f_e_minus2) * len_f_old/len_f_new;

right_side_e_minus2_mu  = interp1(f(f_MI:end),microns(f_MI:end),f_e_minus2);
right_side_e_minus2_pxl = interp1(f(f_MI:end),f_MI:length(f),f_e_minus2)* len_f_old/len_f_new;

e_minus2_mu  =  abs(right_side_e_minus2_mu  - left_side_e_minus2_mu);
e_minus2_pxl =  abs(right_side_e_minus2_pxl - left_side_e_minus2_pxl);



min_i_mu = find(f==min(f));  min_i_mu=min_i_mu(1);
max_i_mu = find(f==max(f));  max_i_mu=max_i_mu(1);

half_f = f(f_MI:end);

diff_f = half_f(2:end) - half_f(1:end-1);
change_Direction = half_f(diff_f>0);
% [~,indx] = min((half_f(2:end))-(half_f(1:end-1)));
First_z_mu =  f_MI + find(half_f == change_Direction(1));
First_z_pxl =  (f_MI + find(half_f == change_Direction(1))) * len_f_old/len_f_new;
% First_z =  f_MI + indx;
zero_x_mu = 2*abs(microns(max_i_mu)-microns(First_z_mu));
zero_x_pxl = 2*abs(pxl(max_i_pxl)-pxl(round(First_z_pxl)));

end

function [theta,R_centers] = Radius_of_convergens(mean_x_center, center_x, mean_y_center, center_y)
X = (center_x-mean_x_center);
Y = (center_y-mean_y_center);
R_centers = sqrt(X.^2+Y.^2);
theta = atan2(Y,X);
end

function [microns] = To_microns (pixels)
global Pixel_size effective_magnification
microns = pixels*Pixel_size/effective_magnification;
end

function totalSum = sum_pixels_in_ellipse(Image, center_point, radii)
%SUM_PIXELS_IN_ELLIPSE Calculates the sum of pixel values within an ellipse.
%   totalSum = sum_pixels_in_ellipse(image, center_point, radii) takes a
%   grayscale image, a center point [center_x, center_y], and ellipse
%   radii [radius_x, radius_y], and returns the sum of pixel values within
%   the ellipse.

[rows, cols] = size(Image);
center_x = center_point(1);
center_y = center_point(2);
radius_x = radii(1);
radius_y = radii(2);
totalSum = 0;
% totalSum_my = 0;
% Iterate through each pixel in the image
for y = 1:rows
    for x = 1:cols
        % Check if the pixel is inside the ellipse
        % (x - center_x)^2 / radius_x^2 + (y - center_y)^2 / radius_y^2 <= 1
        if ((x - center_x)^2 / radius_x^2) + ((y - center_y)^2 / radius_y^2) <= 1
            % totalSum = totalSum + double(Image(y, x));
        end

        if (abs(x - center_x) <= radius_x) && (abs(y - center_y) <= radius_y)
            totalSum = totalSum + double(Image(y, x));
        end

    end
end

end

function [radius_pxl] = Seventy_pers_radius (Image, Full_image,center_point)

[rows, cols] = size(Image);
Full_image_sum = sum(sum(Full_image));
center_x = center_point(1);
center_y = center_point(2);
Ratio = 0;
radius_pxl = 1;

while Ratio<=0.7
    Spot_Sum = 0;
    for y = 1:rows
        for x = 1:cols
            % Check if the pixel is inside the ellipse
            % (x - center_x)^2 / radius_x^2 + (y - center_y)^2 / radius_y^2 <= 1
            if ((x - center_x)^2 / radius_pxl^2) + ((y - center_y)^2 / radius_pxl^2) <= 1
                Spot_Sum = Spot_Sum + double(Image(y, x));
            end
            %
            % if (abs(x - center_x) <= radius_pxl) && (abs(y - center_y) <= radius_pxl)
            %     Spot_Sum = Spot_Sum + double(Image(y, x));
            % end
        end
    end
    Ratio =Spot_Sum/Full_image_sum;
    radius_pxl = radius_pxl+1;
end

end
