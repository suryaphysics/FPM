%%% Code by Surya Kamal to provide proof of concept for Fourier pytchographic microscopy (FPM) %%%

%MIT License

% Copyright (c) 2017 Surya Kamal

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


%%% Initial parameters %%%
inputImage = imread('cameraman.tif'); % Input image
numberOfMasks = 9; % Number of masks
radius = 64; % Size of the mask
center = zeros(numberOfMasks,1); % Array to store the Center of the masks in (x + iy) form
center(1,1) = 0+0i; % Center coordinates of first mask

count = 2;

for theta=0:pi/4:7*(pi/4)   
    center(count,1) = radius*exp(1i*theta); % Generating the coordinates of center of the masks
    count = count + 1;  
end

[x,y] = size(inputImage);
mask = zeros(x,y,numberOfMasks); % 3-D array to store all the masks
lowRes = zeros(x,y,numberOfMasks);
[X,Y] = meshgrid(1:x,1:y);

for i=1:numberOfMasks
    mask(:,:,i) = sqrt((X-(x/2+real(center(i,1)))).^2 + (Y-(y/2+imag(center(i,1)))).^2) < radius; % Storing masks
    imshow(mask(:,:,i),[]);
    lowRes(:,:,i)= abs(fftshift(ifft2(ifftshift((fftshift(fft2(ifftshift(inputImage)))).*mask(:,:,i))))); % Generating LR
    imshow(lowRes(:,:,i),[]);
end


initialSpect = fftshift(fft2(ifftshift(lowRes(:,:,1)))); % broadSpectrum, can begin with any image.
initialNew = initialSpect;
c=1;
% FPM Algorithm
while c<2  % Repeating the whole process 2 times
    for i = 1:numberOfMasks
        
        masked = initialNew.*mask(:,:,i);
        inverse = (fftshift(ifft2(ifftshift(masked))));
        inverseNew = abs(lowRes(:,:,i)) .* exp(1i*angle(inverse));
        newSpectrum = fftshift(fft2(ifftshift(inverseNew)));
        
        for a=1:x
            for b=1:y
                if mask(a,b,i)~=0
                    initialNew(a,b) = newSpectrum(a,b); % Updating 
                end
            end
        end
    end
    c=c+1;
    final = abs(fftshift(ifft2(ifftshift(initialNew))));   
end

 figure;imshow(log(abs(initialNew)),[]); title('Final spectrum');
 figure;imshowpair(lowRes(:,:,1),final,'montage'); title('Low Res Input (left)    High Res Output (right)');
 figure;imshow(angle(fftshift(ifft2(ifftshift(initialNew)))));title('Recovered information')
