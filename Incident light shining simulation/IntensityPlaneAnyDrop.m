function IntMap = IntensityPlaneAnyDrop( rayData, imgDim, res )
%This is for the parallelized collection of intensity, with the hope of
%being able to collect angular data, intensity data all in parallel for a
%fully characterized single experiment. 

%But! This is for the internal droplet, and not just an intensity plane, 
%and is thus specialized for this.

%In multipartrace, we grab the rays final data, location, direction,
%and amplitude. These are then used in order to place them over a matrix
%below with improfile, which then rounds them in order to estimate their
%path and thus 'intensity'. 

%The main goal is to reduce either memory writing of ray cells for large
%ray simulations, and also to reduce expensive computation of other methods

%rayData is a n by 5 array where (x loc) (y loc) (x loc2) (y loc2) (amplitude)

xs = imgDim(1)/(2);
ys = imgDim(2)/(2);
img = zeros(imgDim);
img2 = img;
ext = imgDim(1);
for i = 1:size(rayData,1)
    if rayData(i,5) > 0
        x1 = rayData(i,1)*res + xs;
        x2 = rayData(i,3)*res + xs;
        y1 = rayData(i,2)*res + ys;
        y2 = rayData(i,4)*res + ys;
        a = abs(rayData(i,5));
		
        [cx,cy,~] = improfile(img,[x1 x2],[y1 y2],'bicubic');
        data = round([cx cy]);
		%Below just gets rid of all data which may have just fell outside of the image
        data = data(data(:,1) > 0, :);
        data = data(data(:,2) > 0, :);
        data = data(data(:,2) <= imgDim(1), :);
        data = data(data(:,1) <= imgDim(2), :);
		
		%For each pixel, add the amplitude
        for k = 1:(size(data,1)) %for k = 1:(size(data,1) -1)
            img2(data(k,2),data(k,1)) = img2(data(k,2),data(k,1)) + a;
        end
    end
end
IntMap = img2;
end

