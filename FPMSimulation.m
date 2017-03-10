


inputImage = imread('cameraman.tif');
numberOfMasks = 9;
radius = 64;
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


initialSpect = fftshift(fft2(ifftshift(lowRes(:,:,9)))); % broadSpectrum
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

 figure;imshow(log(abs(initialNew)),[]);
 figure;imshowpair(lowRes(:,:,1),final,'montage');
 figure;imshow(angle(fftshift(ifft2(ifftshift(initialNew)))));