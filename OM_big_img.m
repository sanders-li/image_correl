%{
Sanders Li
Optical Method - Digital Image Correlation
sandersl@usc.edu
Revision History
Date    Changes               Programmer
-----------------------------------------
4/10    Created Program       Sanders Li
%}

clear; close all; clc;
animate = true;
createMovie = false;
setup_mode = 0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  specify the working directory where the a/b image pair is stored
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
working_dir = 'F:\images\Compressed\Rotation1\';  % keep the "\" at the end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the filenames of the a/b image pair in working_dir.  This will create
% a structure where file names are images(1).name for img_a and,  
% images(2).name for img_b.
images = dir([working_dir,'img*.jpg']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Load img_a and img_b and Pre-Process image pair:
%  flatten to 1D (grayscale) if the images are multidimensional (e.g., RGB)
%  This block uses the function 'flatten_rgb_image' which we provided; make
%  sure this function is in your Matlab path.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if createMovie == true
    movie = VideoWriter('visualizer.avi');
    open(movie);
end;

for i = 1:length(images)
    % load image data into a temporary variable
    temp = imread([working_dir,images(i).name]);
    % show original image
    if setup_mode == 1
        figure(1);  subplot(1,2,i);
        imshow(temp);  title(images(i).name,'Interpreter','none'); hold on;
    end

	% Assign flattened (grayscale) image data to the structure "IMAGES"; 
    % this data will be used for interrogation.
    images(i).data = flatten_rgb_image(temp);
end
% End of image Pre-Processing setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Correlation Parameters for the image pair:
%  image box and search box dimensions (in pixels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[img_m,img_n] = size(images(2).data);

Bx = floor(img_m/8);
By = floor(img_n/8); %size of A and B

Sx = floor(Bx*1.75);
Sy = floor(By*1.75); %search box size

startx = floor(Sx/2);
starty = floor(Sy/2); %where to start A, min 1 (start at index 1)

shiftx = Sx;
shifty = Sy;

nx = floor((img_m-Sx-startx)/shiftx)+1; %number of times to shift x
ny = floor((img_n-Sy-starty)/shifty)+1; %number of times to shift y

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  For each image box in img_a, calculate the displacement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nx*ny = total number of displacement calculations (grid points).  This is
% a function of the image size and the Correlation Parameters from above.

hFig = figure(2);
set(hFig, 'Position', [500 200 1280 720]);
hold on;
subplot(3,3,1); imshow(uint8(images(1).data)); 
title(['All A in ' images(1).name],'Interpreter','none'); hold on;
subplot(3,3,2); imshow(uint8(images(1).data)); 
title(['Current A in ' images(1).name],'Interpreter','none'); hold on;
subplot(3,3,4); imshow(uint8(images(2).data)); 
title(['All S in ' images(2).name],'Interpreter','none'); hold on;
subplot(3,3,5); imshow(uint8(images(2).data)); 
title(['Current S and B in ' images(2).name],'Interpreter','none'); hold on;
subplot(3,3,7); imshow(uint8(images(2).data)); 
title('Vectors','Interpreter','none'); hold on;

bdraw = [];
Abox = [];
Sbox = [];

x = zeros(nx,ny);
y = zeros(nx,ny);
u = zeros(nx,ny);
v = zeros(nx,ny);

for p = 1:nx %y-dir for S
    % progress indicator...
    clc;
    fprintf('Progress: %.0f%%... \n',100*(p-1)/nx) %change from nx to ny
    for q = 1:ny %x-dir for S
        fprintf('\t%.0f%%... for row\n',100*(q-1)/ny)
        
        C = zeros(Sx-Bx+1,Sy-By+1);
        C = C-1; %create array of -1s
        
        % pixel array A
        Axpos = startx+Bx/4+(p-1)*shiftx;
        Aypos = starty/2+By/4+(q-1)*shifty;
        A = double(images(1).data(Axpos:Axpos+Bx,Aypos:Aypos+By));  
        % specify array indices and convert to a double
        % NOTE: imshow does not like doubles, 
        % so imshow(uint8(A)) will display A nicely
        A_avg = 1/(Bx*By)*sum(sum(A));  
       
        figure(2); subplot(3,3,3);  imshow(uint8(A)); title('A');
        figure(2); subplot(3,3,1); 
        rectangle('Position',[Aypos-1,Axpos-1,By,Bx],...
            'EdgeColor','r'); 
        %{ 
        %show A in S
        figure(2); subplot(3,3,4); 
        rectangle('Position',[Aypos-1,Axpos-1,By,Bx],...
            'EdgeColor','r'); 
        %}
        figure(2); subplot(3,3,2); delete(Abox);
        Abox = rectangle('Position',[Aypos-1,Axpos-1,By,Bx],...
            'EdgeColor','r'); 
        if animate == true
            drawnow;
        end;
        
        % Find the displacement of A by correlating this pixel array with all 
        % possible destinations B(K,L) in search box S of img_b.
        Sxpos = startx/2+(p-1)*shiftx;
        Sypos = starty/2+(q-1)*shifty;
        S = double(images(2).data(Sxpos:Sxpos+Sx,Sypos:Sypos+Sy));
        
        figure(2); subplot(3,3,4); 
        rectangle('Position',[Sypos-1,Sxpos-1,Sy,Sx],...
            'EdgeColor','b');
        figure(2); subplot(3,3,5); delete(Sbox);
        Sbox = rectangle('Position',[Sypos-1,Sxpos-1,Sy,Sx],...
            'EdgeColor','b');
        if animate == true
            drawnow;
        end;

        for i = 1:(Sx-Bx+1) % x pixel shift within S
            for j = 1:(Sy-By+1) % y pixel shift within S
                %tic %timer function, used to estimate time
                % pixel array B
                % specify array indices within S and convert to a double
                B = double(S(i:Bx+i,j:By+j)); 
                Bxpos = Sxpos+i-1;
                Bypos = Sypos+j-1;
                delete(bdraw);
                figure(2); subplot(3,3,5); 
                bdraw = rectangle('Position',[Bypos-1,Bxpos-1,By,Bx],...
                    'EdgeColor','r');
                if animate == true
                    drawnow;
                end;
                figure(2); subplot(3,3,6); imshow(uint8(B)); title('B'); 
                if animate == true
                    drawnow;
                end;
                if createMovie == true
                    frame = getframe(gcf);
                    writeVideo(movie,frame);
                end;
                
                % Calculate the correlation coefficient, C, for this pixel array.
                % Evaluate C at all possible locations (index shifts I,J).
                % The best correlation determines the displacement of A into img_b.
				%  Note: Double sum below effectively implements Double Riemann sum across k and l in lecture

                B_avg = 1/(Bx*By)*sum(sum(B));
                C(i, j) = sum(sum( (A - A_avg).*(B - B_avg) ))/...
                          sqrt(sum(sum( (A - A_avg).^2 ))*sum(sum( (B - B_avg).^2 )));
                figure(2); subplot(3,3,8); mesh(C); axis([0 Sy-By+1 0 Sx-Bx+1 -1 1]); title('Correlation','Interpreter','none');
                %time = toc; %extension of the timer function
                %speed = time*(Sx-Bx+1)*(Sy-By+1)*nx*ny/3600
            end % j
        end % i
        
        [maxCval1 maxCrow] = max(C);
        [maxCval2 maxCcol] = max(maxCval1);
        maxCx = maxCrow(maxCcol);
        maxCy = maxCcol;
        %maxCval2
        x(p,q) = (Axpos-1 + (Axpos+Bx-1))/2;
        y(p,q) = (Aypos-1 + (Aypos+By-1))/2;
        u(p,q) = (Sxpos+maxCx-1)-Axpos; %x-displacement
        v(p,q) = (Sypos+maxCy-1)-Aypos; %y-displacement
        
        figure(2); subplot(3,3,7); quiver(y,x,v,u,'Marker','.','Color','r'); drawnow;
 %}
    end % q
end % p
fprintf('Processing complete!\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot vectors on shifted image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
imshow(uint8(images(2).data)); % show image 
hold on;         % hold figure
quiver(y,x,v,u,'Marker','.','Color','r','Autoscale','off');   % plot displacement vectors
if createMovie == true
    close(movie);
end;
figure(1); saveas(gcf,'quiver_plot.png');

xlswrite('results.xlsx',x,'x');
xlswrite('results.xlsx',y,'y');
xlswrite('results.xlsx',u,'u');
xlswrite('results.xlsx',v,'v');